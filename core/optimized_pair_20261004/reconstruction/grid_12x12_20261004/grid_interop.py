"""Owned grid-12x12 resident bridge using only the public fixed-byte C ABI.

Import does not load a native library or create a mathematical object. Each
instance owns one query; only the checked CDLL identity is cached per process.
Wrapper clocks cover active methods, excluding idle time between calls. Native
clocks are nested inside create/run/export/destroy, not additional wall time.
"""
import copy
import ctypes as C
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
_LIBRARY = None
_CONTRACT = None
_MAX = 2**63-1
META = tuple('abi state status N L d P e audit expected_products steps_done hinv_created'.split())
STEP_META = tuple('kind phase digit terminal rows inner cols branch native_products mod_g_calls ks_max_bits ks_max_span_bits'.split())
RESULT_META = tuple('kind index rows cols parts degree bits fields bytes status'.split())
COUNTS = tuple('''creates runs exports steps full_products terminal_products zero_products real_products
one_imag_products gaussian_products ks_calls mod_g_calls compose_schedules vec_calls compose_components
compose_arguments constant_components inflate_components inflate_arguments argument_squares hinv_builds hinv_checks
modular_copy_coefficients modular_to_integer_coefficients wire_input_fields wire_output_fields wire_input_bytes
wire_output_bytes debug_bytes poly_inits poly_clears poly_mat_inits poly_mat_clears ctx_inits ctx_clears
bridge_allocs bridge_frees max_ks_bits max_ks_span_bits max_ks_limbs max_vec_slots'''.split())
NATIVE_SECONDS = tuple('''control parse context_setup owner_copy hinv_setup compose_prepare compose_kernel
inflate jump ks_prepare ks_kernel combine_reduce debug export cleanup total'''.split())
WRAPPER_SECONDS = ('load', 'encode', 'create', 'run', 'export', 'decode', 'destroy', 'result_free', 'total')
STATUS = ('OK', 'SHAPE', 'CAPACITY', 'CANONICAL', 'STATE', 'ALLOC', 'ALGEBRA', 'ABI', 'NOT_AVAILABLE', 'COUNTER')
TRACE, INITIAL_A, LEFT, RIGHT, PRODUCT, JUMP_BEFORE, JUMP_AFTER, FINAL_MATRIX = range(8)


class ResidentInputError(ValueError):
    def __init__(self, reason, message):
        self.reason = reason
        super().__init__(message)


class NativeError(RuntimeError):
    def __init__(self, status, metadata=None):
        self.status, self.metadata = int(status), metadata
        self.label = STATUS[self.status] if 0 <= self.status < len(STATUS) else 'UNKNOWN'
        super().__init__('grid-12x12 native '+self.label+' ('+str(self.status)+')')


def checked_input(N, L, d, P, audit=False):
    """Check the finite scalar envelope before buffers, shifting or native use.

    This is a shape check, not a certificate of the Teichmuller modulus. The
    public query pays the original field/parameter construction separately.
    """
    for name, value, minimum, maximum in (('N', N, 1, 128), ('L', L, 2, 128), ('P', P, 1, 8201)):
        if type(value) is not int or not minimum <= value <= maximum:
            raise ResidentInputError(name, name+' is outside the registered exact-integer range')
    if type(d) is not int or not 1 <= d <= N or N % d:
        raise ResidentInputError('d', 'd must be an exact positive divisor of N')
    if type(audit) is not bool:
        raise ResidentInputError('audit', 'audit must be a bool')
    H = (P+7)//8
    e = N//d
    S = d.bit_length()+d.bit_count()+e.bit_length()+e.bit_count()-4
    g_fields, A_fields = d+1, 2*L*L*d
    debug_bytes = (6*S*L*L*d+2*S*d+2*L*L*d)*H if audit else 0
    if any(v > _MAX for v in (g_fields*H, A_fields*H, debug_bytes)):
        raise ResidentInputError('capacity', 'wire size exceeds the registered LP64 capacity')
    return dict(N=N, L=L, d=d, P=P, e=e, H=H, audit=audit, expected_products=S,
                g_fields=g_fields, A_fields=A_fields, g_bytes=g_fields*H, A_bytes=A_fields*H)


def _contract():
    """Share the stdlib-only inventory definition without a generic build import."""
    global _CONTRACT
    if _CONTRACT is None:
        name = 'd_grid_build_defs'
        module = sys.modules.get(name)
        if module is None:
            spec = importlib.util.spec_from_file_location(name, HERE/'build.py')
            if spec is None or spec.loader is None:
                raise RuntimeError('cannot load the grid-12x12 identity definitions')
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            try:
                spec.loader.exec_module(module)
            except BaseException:
                sys.modules.pop(name, None)
                raise
        if Path(module.__file__).resolve() != HERE/'build.py':
            raise RuntimeError('grid-12x12 inventory module has the wrong source path')
        _CONTRACT = module
    return _CONTRACT


def _library():
    global _LIBRARY
    import flint
    if (flint.__version__, flint.__FLINT_VERSION__, flint.ctx.threads) != ('0.9.0', '3.6.0', 1):
        raise RuntimeError('grid-12x12 requires the reviewed single-thread FLINT installation')
    if (sys.maxsize != _MAX or C.sizeof(C.c_int) != 4 or
            any(C.sizeof(t) != 8 for t in (C.c_long, C.c_ulong, C.c_size_t, C.c_void_p))):
        raise RuntimeError('grid-12x12 requires the reviewed LP64 ABI')
    if _LIBRARY is not None:
        return _LIBRARY
    contract = _contract()
    attempt, receipt_path, artifact = contract.selected_build()
    receipt = json.loads(receipt_path.read_text())
    if (not receipt.get('passed') or receipt.get('status') != 'BUILD_ACCEPTABLE' or
            receipt.get('attempt') != attempt or receipt.get('mathematical_operation_executed') is not False):
        raise RuntimeError('grid-12x12 selected build has no acceptable metadata-only receipt')
    fingerprints = {}

    def fingerprint(path):
        path = Path(path).resolve()
        if path not in fingerprints:
            fingerprints[path] = contract.digest(path)
        return fingerprints[path]

    for category, names in (('scientific_sources_sha256', contract.SCIENTIFIC), ('documents_sha256', contract.DOCUMENTS)):
        if set(receipt[category]) != set(names):
            raise RuntimeError('grid-12x12 selected build has an incomplete source inventory')
        for name in names:
            if fingerprint(HERE/name) != receipt[category][name]:
                raise RuntimeError('changed grid-12x12 source/contract: '+name)
    if (receipt['code_review'] != str(contract.review_path(attempt)) or
            fingerprint(contract.review_path(attempt)) != receipt['code_review_sha256']):
        raise RuntimeError('changed grid-12x12 source review')
    for key, expected in (('premises_sha256', contract.PREMISES), ('primary_sources_sha256', contract.PRIMARY)):
        if receipt[key] != {str(p): v for p, v in expected.items()}:
            raise RuntimeError('grid-12x12 receipt does not bind the exact premises/primary inventory')
        for path, value in expected.items():
            if fingerprint(path) != value:
                raise RuntimeError('changed grid-12x12 premise/primary: '+str(path))
    if receipt['evidence_sha256'] != {str(p): v for p, v in contract.EVIDENCE.items()}:
        raise RuntimeError('grid-12x12 receipt does not bind its inherited evidence')
    # Full inherited artifact evidence is checked by the parent/worker manifest.
    # Do not read the inherited rank audit receipt in this per-process bridge load.
    if (fingerprint(artifact) != receipt['artifact_sha256'] or
            receipt['library_sha256'] != contract.LIBRARY_SHA or
            fingerprint(receipt['linked_library']) != contract.LIBRARY_SHA or
            receipt['headers_receipt_sha256'] != contract.HEADERS_SHA or
            fingerprint(receipt['headers_receipt']) != contract.HEADERS_SHA):
        raise RuntimeError('changed grid-12x12 native artifact/library/header identity')
    headers = json.loads(Path(receipt['headers_receipt']).read_text())['headers_sha256']
    if len(headers) != 177:
        raise RuntimeError('grid-12x12 matching header inventory has the wrong size')
    for name, expected in headers.items():
        if fingerprint(contract.REFERENCE/name) != expected:
            raise RuntimeError('changed grid-12x12 matching header: '+name)
    for name, expected in receipt['dependencies_sha256'].items():
        if fingerprint(name) != expected:
            raise RuntimeError('changed grid-12x12 native dependency: '+name)
    lib = C.CDLL(str(artifact), mode=os.RTLD_LOCAL | os.RTLD_NOW)
    for name in ('header_version', 'runtime_version'):
        function = getattr(lib, 'qcb_resident_'+name)
        function.argtypes, function.restype = [], C.c_char_p
    lib.qcb_resident_abi.argtypes, lib.qcb_resident_abi.restype = [], C.c_int
    lib.qcb_resident_slong_size.argtypes, lib.qcb_resident_slong_size.restype = [], C.c_size_t
    if (lib.qcb_resident_abi(), lib.qcb_resident_header_version(), lib.qcb_resident_runtime_version(),
            lib.qcb_resident_slong_size()) != (1, b'3.6.0', b'3.6.0', 8):
        raise RuntimeError('grid-12x12 native constant ABI mismatch')
    pointer = C.c_void_p
    out = C.POINTER(pointer)
    lib.qcb_resident_create.argtypes = [C.c_long, C.c_long, C.c_long, C.c_long, C.c_int,
        pointer, C.c_size_t, pointer, C.c_size_t, out]
    lib.qcb_resident_create.restype = C.c_int
    lib.qcb_resident_run.argtypes, lib.qcb_resident_run.restype = [pointer], C.c_int
    lib.qcb_resident_export.argtypes, lib.qcb_resident_export.restype = [pointer, C.c_int, C.c_long, out], C.c_int
    lib.qcb_resident_destroy.argtypes, lib.qcb_resident_destroy.restype = [out, out], C.c_int
    for prefix in ('', 'result_'):
        for key, result_type in (('meta', C.c_long), ('count', C.c_long), ('seconds', C.c_double)):
            function = getattr(lib, 'qcb_resident_'+prefix+key)
            function.argtypes, function.restype = [pointer, C.c_char_p], result_type
    lib.qcb_resident_step_meta.argtypes, lib.qcb_resident_step_meta.restype = [pointer, C.c_long, C.c_char_p], C.c_long
    lib.qcb_resident_result_data.argtypes, lib.qcb_resident_result_data.restype = [pointer], pointer
    lib.qcb_resident_result_size.argtypes, lib.qcb_resident_result_size.restype = [pointer], C.c_size_t
    lib.qcb_resident_result_free.argtypes, lib.qcb_resident_result_free.restype = [pointer], None
    _LIBRARY = lib
    return lib


def _wire(N, L, d, P, g, A, audit):
    layout = checked_input(N, L, d, P, audit)
    if type(g) is not list or len(g) != d+1:
        raise ResidentInputError('g_shape', 'g must contain exactly d+1 coefficients')
    if any(type(v) is not int or v < 0 or v.bit_length() > P for v in g):
        raise ResidentInputError('g_canonical', 'g coefficients must be canonical Python integers')
    if g[-1] != 1:
        raise ResidentInputError('g_monic', 'g must be monic')
    if type(A) is not list or len(A) != L:
        raise ResidentInputError('A_shape', 'A must have the fixed L by L Gaussian layout')
    for row in A:
        if type(row) is not list or len(row) != L:
            raise ResidentInputError('A_shape', 'A must have the fixed L by L Gaussian layout')
        for pair in row:
            if type(pair) is not list or len(pair) != 2:
                raise ResidentInputError('A_shape', 'each A entry must contain two components')
            for polynomial in pair:
                if type(polynomial) is not list or len(polynomial) != d:
                    raise ResidentInputError('A_shape', 'each A component must contain d fixed coefficients')
                if any(type(v) is not int or v < 0 or v.bit_length() > P for v in polynomial):
                    raise ResidentInputError('A_canonical', 'A coefficients must be canonical Python integers')
    H = layout['H']
    g_wire = b''.join(v.to_bytes(H, 'little') for v in g)
    A_wire = b''.join(v.to_bytes(H, 'little') for row in A for pair in row for poly in pair for v in poly)
    return layout, g_wire, A_wire


def _metadata(lib, handle, receipt=False):
    prefix = 'qcb_resident_result_' if receipt else 'qcb_resident_'
    mg, cg, sg = (getattr(lib, prefix+key) for key in ('meta', 'count', 'seconds'))
    meta = {k: int(mg(handle, k.encode('ascii'))) for k in META}
    counts = {k: int(cg(handle, k.encode('ascii'))) for k in COUNTS}
    seconds = {k: float(sg(handle, k.encode('ascii'))) for k in NATIVE_SECONDS}
    if any(v < 0 for v in (*meta.values(), *counts.values())):
        raise ArithmeticError('grid-12x12 native metadata contains an unknown or negative value')
    if any(not math.isfinite(v) or v < 0 for v in seconds.values()):
        raise ArithmeticError('grid-12x12 native clocks are invalid')
    if meta['abi'] != 1 or meta['state'] not in (1, 2, 3, 4) or meta['status'] >= len(STATUS):
        raise ArithmeticError('grid-12x12 native state/ABI metadata is invalid')
    layout = checked_input(meta['N'], meta['L'], meta['d'], meta['P'], bool(meta['audit']))
    if (meta['audit'] not in (0, 1) or meta['e'] != layout['e'] or
            meta['expected_products'] != layout['expected_products'] or
            meta['steps_done'] > meta['expected_products']):
        raise ArithmeticError('grid-12x12 native schedule metadata contradicts its input')
    return dict(meta=meta, counts=counts, seconds=seconds)


class ResidentReadout:
    """One owned native query. Always close it, preferably with a with block."""

    def __init__(self, N, L, d, P, g, A, audit=False):
        started = perf_counter()
        self._owner = C.c_void_p()
        self._lib, self._receipt = None, None
        self._steps = []
        self._wrapper = dict(seconds={k: 0.0 for k in WRAPPER_SECONDS},
                             result_transfers=0, result_free_calls=0)
        try:
            t = perf_counter()
            self.layout, g_wire, A_wire = _wire(N, L, d, P, g, A, audit)
            # ctypes keeps these real allocated buffers alive for the C call.
            g_buffer, A_buffer = C.create_string_buffer(g_wire), C.create_string_buffer(A_wire)
            self._wrapper['seconds']['encode'] += perf_counter()-t
            t = perf_counter()
            self._lib = _library()
            self._wrapper['seconds']['load'] += perf_counter()-t
            t = perf_counter()
            status = int(self._lib.qcb_resident_create(N, L, d, P, int(audit),
                g_buffer, len(g_wire), A_buffer, len(A_wire), C.byref(self._owner)))
            self._wrapper['seconds']['create'] += perf_counter()-t
            if status:
                if self._owner.value:
                    self.close()
                    raise ArithmeticError('grid-12x12 failed create returned a live owner')
                raise NativeError(status)
            if not self._owner.value:
                raise ArithmeticError('grid-12x12 successful create returned no owner')
        finally:
            self._wrapper['seconds']['total'] += perf_counter()-started

    def _live(self):
        if not self._owner.value:
            raise ResidentInputError('closed', 'grid-12x12 owner is closed')

    def _owner_snapshot(self):
        self._live()
        snapshot = _metadata(self._lib, self._owner)
        steps = []
        for i in range(snapshot['meta']['steps_done']):
            row = {k: int(self._lib.qcb_resident_step_meta(self._owner, i, k.encode('ascii'))) for k in STEP_META}
            if any(v < 0 for v in row.values()):
                raise ArithmeticError('grid-12x12 step metadata contains an unknown value')
            steps.append(row)
        snapshot['steps'] = steps
        return snapshot

    def run(self):
        started = perf_counter()
        try:
            self._live()
            t = perf_counter()
            status = int(self._lib.qcb_resident_run(self._owner))
            self._wrapper['seconds']['run'] += perf_counter()-t
            if status:
                raise NativeError(status, self._owner_snapshot())
            return self
        finally:
            self._wrapper['seconds']['total'] += perf_counter()-started

    def _free_result(self, result):
        if result.value:
            t = perf_counter()
            self._lib.qcb_resident_result_free(result)
            result.value = None
            self._wrapper['result_free_calls'] += 1
            self._wrapper['seconds']['result_free'] += perf_counter()-t

    def export(self, kind, index=0):
        started = perf_counter()
        result = C.c_void_p()
        try:
            self._live()
            if type(kind) is not int or not -(2**31) <= kind < 2**31:
                raise ResidentInputError('kind', 'export kind must fit a C int')
            if type(index) is not int or not -_MAX-1 <= index <= _MAX:
                raise ResidentInputError('index', 'export index must fit a C long')
            t = perf_counter()
            status = int(self._lib.qcb_resident_export(self._owner, kind, index, C.byref(result)))
            self._wrapper['seconds']['export'] += perf_counter()-t
            if result.value:
                self._wrapper['result_transfers'] += 1
            if status:
                if result.value:
                    raise ArithmeticError('grid-12x12 rejected export returned a result')
                raise NativeError(status, self._owner_snapshot())
            if not result.value:
                raise ArithmeticError('grid-12x12 successful export returned no result')
            t = perf_counter()
            meta = {k: int(self._lib.qcb_resident_result_meta(result, k.encode('ascii'))) for k in RESULT_META}
            length = int(self._lib.qcb_resident_result_size(result))
            address = self._lib.qcb_resident_result_data(result)
            L, d, P, H = (self.layout[k] for k in ('L', 'd', 'P', 'H'))
            if (meta['kind'] != kind or meta['index'] != index or meta['status'] != 0 or
                    meta['degree'] != d or meta['bits'] != P or
                    meta['rows'] not in (1, L) or meta['cols'] not in (1, L) or
                    meta['parts'] not in (1, 2)):
                raise ArithmeticError('grid-12x12 export metadata contradicts its fixed layout')
            if kind == TRACE:
                shape = (1, 1, 2)
            elif kind in (INITIAL_A, LEFT, RIGHT, FINAL_MATRIX):
                shape = (L, L, 2)
            elif kind in (JUMP_BEFORE, JUMP_AFTER):
                shape = (1, 1, 1)
            elif kind == PRODUCT:
                terminal = int(self._lib.qcb_resident_step_meta(self._owner, index, b'terminal'))
                if terminal not in (0, 1):
                    raise ArithmeticError('grid-12x12 product has no valid terminal metadata')
                shape = (1, 1, 2) if terminal else (L, L, 2)
            else:
                raise ArithmeticError('grid-12x12 native accepted an unknown export kind')
            if (meta['rows'], meta['cols'], meta['parts']) != shape:
                raise ArithmeticError('grid-12x12 export has the wrong complete mathematical shape')
            fields = meta['rows']*meta['cols']*meta['parts']*d
            if meta['fields'] != fields or meta['bytes'] != fields*H or length != fields*H or not address:
                raise ArithmeticError('grid-12x12 export byte count contradicts its layout')
            wire = C.string_at(address, length)
            position = 0
            coefficients = []
            for _ in range(meta['rows']):
                row = []
                for _ in range(meta['cols']):
                    pair = []
                    for _ in range(meta['parts']):
                        polynomial = []
                        for _ in range(d):
                            value = int.from_bytes(wire[position:position+H], 'little')
                            position += H
                            if value.bit_length() > P:
                                raise ArithmeticError('grid-12x12 export contains a noncanonical coefficient')
                            polynomial.append(value)
                        pair.append(polynomial)
                    row.append(pair)
                coefficients.append(row)
            self._wrapper['seconds']['decode'] += perf_counter()-t
            return dict(meta=meta, coefficients=coefficients)
        finally:
            self._free_result(result)
            self._wrapper['seconds']['total'] += perf_counter()-started

    def close(self):
        """Destroy once; retain an independent Python receipt, including cleanup."""
        if not self._owner.value:
            return self._receipt
        started = perf_counter()
        result = C.c_void_p()
        try:
            # The final C receipt has owner counters but no separate step getter.
            # Copy the committed finite step list before destroying its owner.
            before, snapshot_error = None, None
            try:
                before = self._owner_snapshot()
                self._steps = before['steps']
            except Exception as error:
                # Even malformed native diagnostics must not bypass destroy.
                snapshot_error = error
            t = perf_counter()
            status = int(self._lib.qcb_resident_destroy(C.byref(self._owner), C.byref(result)))
            self._wrapper['seconds']['destroy'] += perf_counter()-t
            if result.value:
                self._wrapper['result_transfers'] += 1
            if status:
                raise NativeError(status, before)
            if self._owner.value or not result.value:
                raise ArithmeticError('grid-12x12 destroy did not transfer its final receipt')
            if snapshot_error is not None:
                raise snapshot_error
            t = perf_counter()
            if (self._lib.qcb_resident_result_meta(result, b'kind') != -1 or
                    self._lib.qcb_resident_result_meta(result, b'fields') != 0 or
                    self._lib.qcb_resident_result_meta(result, b'bytes') != 0 or
                    self._lib.qcb_resident_result_size(result) != 0 or
                    self._lib.qcb_resident_result_data(result)):
                raise ArithmeticError('grid-12x12 destroy receipt unexpectedly has a payload')
            final = _metadata(self._lib, result, receipt=True)
            final['steps'] = self._steps
            if (final['meta']['state'], final['meta']['status'], final['meta']['steps_done']) != (
                    before['meta']['state'], before['meta']['status'], before['meta']['steps_done']):
                raise ArithmeticError('grid-12x12 destroy changed the completed mathematical state')
            self._receipt = final
            self._wrapper['seconds']['decode'] += perf_counter()-t
            return self._receipt
        finally:
            self._free_result(result)
            self._wrapper['seconds']['total'] += perf_counter()-started

    def statistics_snapshot(self):
        native = self._owner_snapshot() if self._owner.value else self._receipt
        if native is None:
            raise ResidentInputError('closed', 'grid-12x12 has no successful owner or final receipt')
        # This observation is pure: repeated snapshots after close are equal.
        # Getter/copy cost is paid by the enclosing query/validation wall, not
        # retroactively added to the already finalized wrapper receipt.
        return dict(native=copy.deepcopy(native), wrapper=copy.deepcopy(self._wrapper))

    def __enter__(self):
        self._live()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False
