"""D3 FLINT-only batch conversion. No build, load or math on import."""
import ctypes as C
import hashlib
import json
import os
from pathlib import Path
from time import perf_counter

from flint import fmpz

HERE = Path(__file__).resolve().parent
C_SOURCE_ROOT = HERE.parent
_LIBRARY = None
_SOURCES = ('compose.c', 'interop.py', 'build.py')


def _digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _library():
    global _LIBRARY
    import flint
    if flint.ctx.threads != 1:
        raise RuntimeError('registered batch-composition requires one thread')
    if _LIBRARY is not None:
        return _LIBRARY
    if (flint.__version__, flint.__FLINT_VERSION__) != ('0.9.0', '3.6.0'):
        raise RuntimeError('unreviewed batch-composition version')
    receipt = json.loads((C_SOURCE_ROOT/'build/receipt.json').read_text())
    artifact = C_SOURCE_ROOT/'build/libqcb_batch.so'
    if not receipt['passed']:
        raise RuntimeError('unaccepted batch-composition build')
    for name in _SOURCES:
        if receipt['sources_sha256'][name] != _digest(C_SOURCE_ROOT/name):
            raise RuntimeError('stale batch-composition source: '+name)
    if receipt['protocol_sha256'] != _digest(C_SOURCE_ROOT/'PROTOCOL.md'):
        raise RuntimeError('batch-composition protocol changed')
    if receipt['artifact_sha256'] != _digest(artifact):
        raise RuntimeError('batch-composition artifact changed')
    if receipt['library_sha256'] != _digest(Path(receipt['linked_library'])):
        raise RuntimeError('bundled FLINT library changed')
    if receipt['headers_receipt_sha256'] != _digest(Path(receipt['headers_receipt'])):
        raise RuntimeError('reviewed header receipt changed')
    lib = C.CDLL(str(artifact), mode=os.RTLD_LOCAL | os.RTLD_NOW)
    for name in ('header_version', 'runtime_version'):
        getattr(lib, 'qcb_batch_'+name).restype = C.c_char_p
    lib.qcb_batch_slong_size.restype = C.c_size_t
    if (lib.qcb_batch_header_version(), lib.qcb_batch_runtime_version()) != (b'3.6.0', b'3.6.0'):
        raise RuntimeError('batch-composition header/runtime mismatch')
    if lib.qcb_batch_abi() != 1 or lib.qcb_batch_slong_size() != C.sizeof(C.c_long):
        raise RuntimeError('batch-composition ABI mismatch')
    lib.qcb_batch_create.argtypes = [C.c_long, C.c_char_p, C.c_char_p, C.POINTER(C.c_void_p)]
    lib.qcb_batch_create.restype = C.c_int
    lib.qcb_batch_destroy.argtypes = [C.c_void_p]
    lib.qcb_batch_destroy.restype = None
    lib.qcb_batch_setup_seconds.argtypes = [C.c_void_p]
    lib.qcb_batch_setup_seconds.restype = C.POINTER(C.c_double)
    lib.qcb_batch_compose.argtypes = [C.c_void_p, C.c_long, C.POINTER(C.c_char_p), C.c_char_p,
                                     C.POINTER(C.c_void_p)]
    lib.qcb_batch_compose.restype = C.c_int
    for name, kind in (('strings', C.c_void_p), ('dimensions', C.c_long), ('seconds', C.c_double)):
        getter = getattr(lib, 'qcb_batch_result_'+name)
        getter.argtypes = [C.c_void_p]
        getter.restype = C.POINTER(kind)
    lib.qcb_batch_result_free.argtypes = [C.c_void_p]
    lib.qcb_batch_result_free.restype = None
    _LIBRARY = lib
    return lib


def _status(status, operation):
    if status == 5:
        raise MemoryError('batch-composition adapter allocation failed')
    if status:
        description = {1: 'invalid input', 2: 'allocation-size overflow',
                       3: 'invalid polynomial serialization',
                       4: 'preinverse verification failed'}.get(status, 'unexpected native status')
        raise ValueError(operation+': '+description+' ('+str(status)+')')


def _encode(poly):
    coefficients = list(poly)
    if not coefficients:
        return b'0'
    return (str(len(coefficients))+'  '+' '.join(fmpz(int(x)).str(base=10, condense=0) for x in coefficients)).encode('ascii')


def _decode(pointer, ctx):
    parts = C.string_at(pointer).split()
    if not parts:
        raise ArithmeticError('empty native batch polynomial')
    count = int(parts[0])
    if count != len(parts)-1:
        raise ArithmeticError('native batch polynomial length mismatch')
    return ctx([fmpz(x.decode('ascii')) for x in parts[1:]])


class BatchComposer:
    """One query-owned g/inverse; every compose call rebuilds its powers table."""

    def __init__(self, R, *, inverse_override=None):
        self._handle = C.c_void_p()
        self._lib = None
        self.R = R
        self._statistics = dict(
            calls=0, native_calls=0, polynomials=0, hinv_builds=0, hinv_checks=0,
            inverse_overrides=0, zero_argument_calls=0, empty_calls=0,
            max_batch=0, max_s=0, max_k=0, max_native_slots=0,
            temp_result_bits_bound=0, input_bytes=0, output_bytes=0,
            total_seconds=0.0, load_seconds=0.0, setup_seconds=0.0,
            context_parse_seconds=0.0, hinv_seconds=0.0, hinv_check_seconds=0.0,
            input_conversion_seconds=0.0, native_parse_seconds=0.0,
            native_compose_seconds=0.0, native_serialize_seconds=0.0,
            native_cleanup_seconds=0.0, output_conversion_seconds=0.0,
            cleanup_seconds=0.0)
        started = perf_counter()
        self.bits = int(R.mod).bit_length()-1
        if R.d < 2 or self.bits < 1 or R.mod != 1 << self.bits:
            raise ValueError('BatchComposer requires d>=2 and modulus 2^P')
        if int(R.ctx.modulus()) != R.mod:
            raise ValueError('inconsistent output coefficient context')
        if self.bits > ((1 << (8*C.sizeof(C.c_long)-1))-1)//4:
            raise ValueError('precision exceeds native ABI')
        self._check_poly(R.g, R.d+1)
        if R.g.degree() != R.d or int(R.g[R.d]) != 1:
            raise ValueError('BatchComposer requires a monic degree-d modulus')
        if inverse_override is not None:
            self._check_poly(inverse_override, R.d+1)
        load_start = perf_counter()
        self._lib = _library()
        load_seconds = perf_counter()-load_start
        self._statistics['load_seconds'] = load_seconds
        try:
            encoded_g = _encode(R.g)
            encoded_inverse = None if inverse_override is None else _encode(inverse_override)
            status = self._lib.qcb_batch_create(self.bits, encoded_g, encoded_inverse, C.byref(self._handle))
            _status(status, 'batch context creation')
            timings = self._lib.qcb_batch_setup_seconds(self._handle)
            for index, name in enumerate(('context_parse_seconds', 'hinv_seconds', 'hinv_check_seconds')):
                self._statistics[name] = timings[index]
            self._statistics['hinv_builds'] = int(inverse_override is None)
            self._statistics['inverse_overrides'] = int(inverse_override is not None)
            self._statistics['hinv_checks'] = 1
            self._statistics['input_bytes'] = len(encoded_g)+(len(encoded_inverse) if encoded_inverse else 0)
            self._statistics['setup_seconds'] = perf_counter()-started-load_seconds
            self._statistics['total_seconds'] = perf_counter()-started
        except BaseException:
            self.close()
            raise

    def _check_poly(self, polynomial, length_limit):
        if not hasattr(polynomial, 'modulus') or not hasattr(polynomial, 'degree'):
            raise TypeError('expected a public FLINT modular polynomial')
        if int(polynomial.modulus()) != self.R.mod:
            raise ValueError('inconsistent coefficient context')
        if polynomial.degree() >= length_limit:
            raise ValueError('polynomial degree exceeds batch contract')

    def compose_many(self, polynomials, argument):
        if not self._handle:
            raise RuntimeError('BatchComposer is closed')
        started = perf_counter()
        stats = self._statistics
        stats['calls'] += 1
        self._check_poly(argument, self.R.d)
        values = list(polynomials)
        for value in values:
            self._check_poly(value, self.R.d)
        if not values:
            stats['empty_calls'] += 1
            stats['total_seconds'] += perf_counter()-started
            return []
        count = len(values)
        long_max = (1 << (8*C.sizeof(C.c_long)-1))-1
        if count >= long_max or self.R.d > long_max//count:
            raise ValueError('batch dimensions exceed native ABI')
        encoded = [_encode(value) for value in values]
        encoded_argument = _encode(argument)
        inputs = (C.c_char_p*count)(*encoded)
        stats['input_bytes'] += sum(map(len, encoded))+len(encoded_argument)
        stats['input_conversion_seconds'] += perf_counter()-started
        result = C.c_void_p()
        try:
            status = self._lib.qcb_batch_compose(self._handle, count, inputs,
                                                encoded_argument, C.byref(result))
            _status(status, 'batch composition')
            dimensions = self._lib.qcb_batch_result_dimensions(result)
            if dimensions[0] != count:
                raise ArithmeticError('native batch result count mismatch')
            for index, name in ((0, 'max_batch'), (1, 'max_s'), (2, 'max_k'),
                                (3, 'max_native_slots'), (4, 'temp_result_bits_bound')):
                stats[name] = max(stats[name], int(dimensions[index]))
            elapsed = self._lib.qcb_batch_result_seconds(result)
            for index, name in enumerate(('native_parse_seconds', 'native_compose_seconds',
                                          'native_serialize_seconds', 'native_cleanup_seconds')):
                stats[name] += elapsed[index]
            converted = perf_counter()
            pointers = self._lib.qcb_batch_result_strings(result)
            output = []
            for j in range(count):
                stats['output_bytes'] += len(C.string_at(pointers[j]))
                polynomial = _decode(pointers[j], self.R.ctx)
                if polynomial.degree() >= self.R.d:
                    raise ArithmeticError('native batch returned unreduced polynomial')
                output.append(polynomial)
            stats['output_conversion_seconds'] += perf_counter()-converted
            stats['native_calls'] += 1
            stats['polynomials'] += count
            stats['zero_argument_calls'] += int(not argument)
            return output
        finally:
            cleaned = perf_counter()
            if result:
                self._lib.qcb_batch_result_free(result)
            stats['cleanup_seconds'] += perf_counter()-cleaned
            stats['total_seconds'] += perf_counter()-started

    def statistics_snapshot(self):
        return dict(self._statistics)

    def close(self):
        if self._handle and self._lib is not None:
            started = perf_counter()
            self._lib.qcb_batch_destroy(self._handle)
            self._handle = C.c_void_p()
            elapsed = perf_counter()-started
            self._statistics['cleanup_seconds'] += elapsed
            self._statistics['total_seconds'] += elapsed

    def __enter__(self):
        if not self._handle:
            raise RuntimeError('BatchComposer is closed')
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass
