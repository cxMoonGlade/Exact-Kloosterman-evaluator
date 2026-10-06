"""D5 owned first-transform bridge; no native loading or math on import."""
import ctypes as C
import hashlib
import json
import os
from pathlib import Path
from time import perf_counter

from flint import fmpz

HERE = Path(__file__).resolve().parent
_LIBRARY = None
_SOURCES = ('guard.c', 'guard_interop.py', 'build.py')
_API = HERE.parent/'first_transform/API_REVIEW.md'
_PREMISES = {
    str(HERE.parent.parent/'PRECISION_REDUCTION_PROPOSAL.md'):
        '540e1b66b57e50139f36cf96b4f668dccef75e89d658de67deb94b02bca7565d',
    str(HERE.parent.parent/'PRECISION_REDUCTION_REVIEW.md'):
        'f88bd8175fb7f73b71fb487d4127724b7e85b87d5d3538583be41c6d16f32e88',
}
_STAGES = ('setup', 'parse', 'kernel', 'tree', 'weights', 'solve',
           'serialize', 'workspace_cleanup')
_COUNTS = ('Q', 'r', 'rhs_count', 'tree_nodes', 'kernel_steps',
           'kernel_odd_inverses', 'unit_inverses', 'weight_inverses',
           'tree_products', 'rhs_convolutions', 'interpolation_products',
           'interpolation_additions', 'node_values', 'output_slots',
           'tree_W_slots', 'tree_C_slots', 'intermediate_coefficient_bits_bound',
           'kernel_term_products')


def _digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _library():
    global _LIBRARY
    import flint
    if flint.ctx.threads != 1:
        raise RuntimeError('D5 requires one native thread')
    if _LIBRARY is not None:
        return _LIBRARY
    if (flint.__version__, flint.__FLINT_VERSION__) != ('0.9.0', '3.6.0'):
        raise RuntimeError('unreviewed D5 Python-FLINT/native version')
    receipt = json.loads((HERE/'build/receipt.json').read_text())
    artifact = HERE/'build/libqcb_guard.so'
    if not receipt['passed']:
        raise RuntimeError('unaccepted D5 build')
    for name in _SOURCES:
        if receipt['sources_sha256'][name] != _digest(HERE/name):
            raise RuntimeError('stale D5 bridge source: '+name)
    if receipt['protocol_sha256'] != _digest(HERE/'PROTOCOL.md'):
        raise RuntimeError('changed D5 protocol')
    if receipt['api_path'] != str(_API) or receipt['api_sha256'] != _digest(_API):
        raise RuntimeError('changed inherited D4 API identity')
    if receipt.get('premises_sha256') != _PREMISES:
        raise RuntimeError('D5 receipt does not bind the reviewed precision proof')
    for path, expected in _PREMISES.items():
        if _digest(Path(path)) != expected:
            raise RuntimeError('changed D5 mathematical premise: '+path)
    if receipt['artifact_sha256'] != _digest(artifact):
        raise RuntimeError('changed D5 native artifact')
    if receipt['library_sha256'] != _digest(Path(receipt['linked_library'])):
        raise RuntimeError('changed bundled FLINT library')
    if receipt['headers_receipt_sha256'] != _digest(Path(receipt['headers_receipt'])):
        raise RuntimeError('changed reviewed header receipt')
    lib = C.CDLL(str(artifact), mode=os.RTLD_LOCAL | os.RTLD_NOW)
    for name in ('header_version', 'runtime_version'):
        getter = getattr(lib, 'qcb_first_'+name)
        getter.argtypes = []
        getter.restype = C.c_char_p
    lib.qcb_first_abi.argtypes = []
    lib.qcb_first_abi.restype = C.c_int
    lib.qcb_first_slong_size.argtypes = []
    lib.qcb_first_slong_size.restype = C.c_size_t
    if (lib.qcb_first_header_version(), lib.qcb_first_runtime_version()) != (b'3.6.0', b'3.6.0'):
        raise RuntimeError('D5 header/runtime mismatch')
    if lib.qcb_first_abi() != 2 or lib.qcb_first_slong_size() != C.sizeof(C.c_long):
        raise RuntimeError('D5 ABI mismatch')
    lib.qcb_first_solve.argtypes = [C.c_long, C.c_long, C.c_char_p,
                                    C.POINTER(C.c_char_p), C.POINTER(C.c_void_p)]
    lib.qcb_first_solve.restype = C.c_int
    lib.qcb_first_result_free.argtypes = [C.c_void_p]
    lib.qcb_first_result_free.restype = None
    for name, kind in (('strings', C.c_void_p), ('counts', C.c_long), ('seconds', C.c_double)):
        getter = getattr(lib, 'qcb_first_result_'+name)
        getter.argtypes = [C.c_void_p]
        getter.restype = C.POINTER(kind)
    _LIBRARY = lib
    return lib


def _status(status):
    if status == 5:
        raise MemoryError('D5 adapter allocation failed')
    if status in (4, 6):
        message = 'odd unit inverse failed' if status == 4 else 'binary division is not exact'
        raise ArithmeticError('D5 '+message)
    if status:
        message = {1: 'invalid input contract', 2: 'dimension/allocation overflow',
                   3: 'invalid polynomial wire'}.get(status, 'unexpected native status')
        raise ValueError('D5 '+message+' ('+str(status)+')')


def _encode(values, Q):
    if not isinstance(values, (list, tuple)) or len(values) > Q:
        raise ValueError('D5 input must be a sequence with at most Q coefficients')
    converted = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, (int, fmpz)):
            raise TypeError('D5 coefficients must be integers')
        value = int(value)
        if value < 0 or value.bit_length() > Q:
            raise ValueError('D5 coefficient is not a canonical Q-bit residue')
        converted.append(fmpz(value).str(base=10, condense=0))
    if not converted:
        return b'0'
    return (str(len(converted))+'  '+' '.join(converted)).encode('ascii')


def _decode(data, Q):
    parts = data.split()
    if not parts:
        raise ArithmeticError('empty D5 result polynomial')
    count = int(parts[0])
    if count < 0 or count > Q or count != len(parts)-1:
        raise ArithmeticError('D5 result coefficient count mismatch')
    output = [int(fmpz(value.decode('ascii'))) for value in parts[1:]]
    if any(value < 0 or value.bit_length() > Q for value in output):
        raise ArithmeticError('D5 returned a noncanonical coefficient')
    return output + [0]*(Q-count)


def first_transform(phi, H, Q, r):
    """Return r shifted-unit solutions then three H solutions, each Q slots.

    The metadata total includes its mutually exclusive native/Python stages;
    tree slot counts and the explicit-arithmetic bit bound are not RSS samples.
    """
    started = perf_counter()
    if type(Q) is not int or Q < 2 or type(r) is not int or not 1 <= r <= Q:
        raise ValueError('D5 requires Q>=2 and 1<=r<=Q')
    long_max = (1 << (8*C.sizeof(C.c_long)-1))-1
    if Q > long_max//16 or Q*(r+3) > long_max or 2*(Q-1)*(r+3) > long_max:
        raise ValueError('D5 dimensions exceed the native ABI')
    # C H[3] is a pointer: enforce three accessible slots before the call.
    if not isinstance(H, (list, tuple)) or len(H) != 3:
        raise ValueError('D5 requires exactly three H polynomials')
    stages = {}
    mark = perf_counter()
    encoded_phi = _encode(phi, Q)
    encoded_h = [_encode(polynomial, Q) for polynomial in H]
    pointers = (C.c_char_p*3)(*encoded_h)
    stages['encode'] = perf_counter()-mark
    input_bytes = len(encoded_phi)+sum(map(len, encoded_h))
    mark = perf_counter()
    lib = _library()
    stages['load'] = perf_counter()-mark
    result = C.c_void_p()
    output_bytes = 0
    try:
        _status(lib.qcb_first_solve(Q, r, encoded_phi, pointers, C.byref(result)))
        if not result:
            raise ArithmeticError('D5 success returned no result')
        mark = perf_counter()
        native_counts = lib.qcb_first_result_counts(result)
        counts = {name: int(native_counts[index]) for index, name in enumerate(_COUNTS)}
        if (counts['Q'], counts['r'], counts['rhs_count']) != (Q, r, r+3):
            raise ArithmeticError('D5 result dimensions disagree')
        native_seconds = lib.qcb_first_result_seconds(result)
        stages.update({name: native_seconds[index] for index, name in enumerate(_STAGES)})
        strings = lib.qcb_first_result_strings(result)
        output = []
        for index in range(r+3):
            if not strings[index]:
                raise ArithmeticError('D5 result string is missing')
            data = C.string_at(strings[index])
            output_bytes += len(data)
            output.append(_decode(data, Q))
        stages['decode'] = perf_counter()-mark
    finally:
        mark = perf_counter()
        if result:
            lib.qcb_first_result_free(result)
        stages['result_cleanup'] = perf_counter()-mark
    metadata = dict(Q=Q, r=r, rhs_count=r+3, input_bytes=input_bytes,
                    output_bytes=output_bytes, stages=stages, counts=counts,
                    precision=dict(kernel_initial=4*Q, convolution=3*Q, tree=2*Q, output=Q),
                    coefficient_bound_scope='explicit coefficient arithmetic; excludes FLINT packing/workspace',
                    algorithm='finite_first_transform_4Q_3Q_2Q')
    metadata['seconds'] = perf_counter()-started
    return output, metadata
