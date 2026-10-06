"""Owned D13 binary bridge. Import does not load native libraries or do math.

Local timers stop before returning; later Python temporary destruction is
paid by the enclosing acquisition and cold process, not this local snapshot.
"""
import ctypes as C
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
_LIBRARY = None
_MAX = sys.maxsize
_SOURCES = ('block.c', 'block_interop.py', 'build.py', 'block_acquisition.py', 'diagnose.py')
_DOCUMENTS = ('README.md', 'PROTOCOL.md', 'API_CONTRACT.md', 'REVIEW.md')
_PREMISES = {
    ROOT/'BLOCK_ELIMINATION_DESIGN_DRAFT.md': '8b07f2f1a2e7db77dcdc316d209de7d7d57b9c0b5cf922c716b050bdd57c8123',
    ROOT/'BLOCK_ELIMINATION_DESIGN_REVIEW.md': '163e7e6506b102a94ecf77d797cba824c6829ca886561c2461196ed621be5609',
    ROOT/'BLOCK_ELIMINATION_API_NOTE.md': '01088dc06238a2b907aa032ce167cee23950bb0e230f994e82f58ae8d042d513',
    ROOT/'BLOCK_ELIMINATION_MATRIX_API_CONTRACT.md': 'eb5193db7a3f659b7411406d218b71def6f5a95e156e91b403bef389b818a2cc',
    ROOT/'BLOCK_ELIMINATION_SCHEDULE_DRAFT.md': 'c0a206c050b37cbf9352bc3934b6caf07b180543a6411496161c7b91644937bc',
    ROOT/'BLOCK_ELIMINATION_SCHEDULE_REVIEW.md': 'b1b388fbee0d8429c7c7c1cbd953afd1c133d9a60a79a443d50a7a3956a2c057',
}
_META = ('status', 'mode', 'Q', 'L', 'n', 'node', 'k', 'next_bits', 'residual_bits', 'audit', 'field_bytes')
_WORDS = dict(char_bits=8, flint_bits=64, long_bytes=8, slong_bytes=8, ulong_bytes=8, size_t_bytes=8)
_COUNTS = tuple('''hat_count rhs_count hat_slots unit_slots normalized_rhs_slots prefix_slots output_slots
fall_divisions fall_products unit_inversions ratio_products contexts rhs_contexts H_filter_slots gamma_valuation_checks
h0_rhs_products h1_rhs_products second_input_matrix_slots native_factorizations
block_groups constant_groups short_groups successor_groups terminal_groups block_coefficient_slots
block_successor_slots taylor_products short_unit_inverses short_power_products triangular_dot_terms
triangular_diagonal_products diagonal_power_products long_basis_products correction_products long_state_products
block_divisions b_integrality_slots remainder_checked_slots power2_checked_slots block_shifted_slots
matrix_inits matrix_clears matrix_slots_initialized matrix_slots_cleared matrix_coefficient_gets
matrix_coefficient_sets matrix_products matrix_input_entries matrix_output_entries audit_next_slots
input_fields_parsed input_bytes_parsed output_fields_written output_bytes_written array_import_calls array_export_calls
scalar_coefficient_gets scalar_coefficient_sets integer_poly_get_calls integer_poly_get_slots
integer_poly_set_calls integer_poly_set_slots bulk_right_shift_calls bulk_right_shift_input_slots
canonical_checked_fields filtered_checked_slots audit_fields_written fmpz_inits fmpz_clears
poly_inits poly_clears context_inits context_clears workspace_alloc_calls workspace_free_calls
workspace_allocated_bytes workspace_freed_bytes'''.split())
_NATIVE_STAGES = ('parse_check', 'contexts', 'RHS', 'units', 'ratios', 'block_prepare',
    'taylor_matmul', 'short_solve', 'long_reconstruct', 'residual_check_shift',
    'assembly', 'serialize', 'workspace_cleanup', 'total')
_STATUS = ('OK', 'ARGUMENT', 'CAPACITY', 'WIRE_LENGTH', 'NONCANONICAL', 'ALLOCATION',
           'GAMMA', 'FILTER', 'HAT_OR_REMAINDER', 'NONINTEGRAL', 'NONUNIT', 'INTERNAL', 'ABI')


class NativeBlockError(RuntimeError):
    """A native rejection with retained actual diagnostics, if a head existed."""
    def __init__(self, status, metadata):
        self.status, self.metadata = status, metadata
        label = _STATUS[status] if 0 <= status < len(_STATUS) else 'UNKNOWN'
        super().__init__('D13 native '+label+' ('+str(status)+')')


def _digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def _library():
    global _LIBRARY
    import flint
    if (flint.__version__, flint.__FLINT_VERSION__, flint.ctx.threads) != ('0.9.0', '3.6.0', 1):
        raise RuntimeError('D13 requires the reviewed single-thread Python/FLINT installation')
    if any(C.sizeof(t) != 8 for t in (C.c_long, C.c_ulong, C.c_size_t, C.c_void_p)) or _MAX != 2**63-1:
        raise RuntimeError('D13 requires the reviewed LP64 Python ABI')
    if _LIBRARY is not None:
        return _LIBRARY
    # One load checks each file's bytes once even when multiple receipt fields
    # name the same library. This local memo is discarded after this load.
    fingerprints = {}

    def fingerprint(path):
        path = path.resolve()
        if path not in fingerprints:
            fingerprints[path] = _digest(path)
        return fingerprints[path]

    receipt = json.loads((HERE/'build/receipt.json').read_text())
    if not receipt.get('passed') or receipt.get('status') != 'BUILD_ACCEPTABLE':
        raise RuntimeError('unaccepted D13 build')
    for name in _SOURCES:
        if receipt['scientific_sources_sha256'][name] != fingerprint(HERE/name):
            raise RuntimeError('changed D13 scientific source: '+name)
    for name in _DOCUMENTS:
        if receipt['documents_sha256'][name] != fingerprint(HERE/name):
            raise RuntimeError('changed D13 contract: '+name)
    if receipt['code_review_sha256'] != fingerprint(HERE/'CODE_REVIEW.md'):
        raise RuntimeError('changed D13 source review')
    expected_premises = {str(path): value for path, value in _PREMISES.items()}
    if receipt['premises_sha256'] != expected_premises:
        raise RuntimeError('D13 build does not bind the finite premises')
    for path, value in _PREMISES.items():
        if fingerprint(path) != value:
            raise RuntimeError('changed D13 premise: '+str(path))
    artifact = HERE/'build/libqcb_block.so'
    if (receipt['artifact_sha256'] != fingerprint(artifact)
            or receipt['library_sha256'] != fingerprint(Path(receipt['linked_library']))
            or receipt['headers_receipt_sha256'] != fingerprint(Path(receipt['headers_receipt']))):
        raise RuntimeError('changed D13 native build identity')
    for path, expected in receipt['dependencies_sha256'].items():
        if fingerprint(Path(path)) != expected:
            raise RuntimeError('changed D13 native dependency: '+path)
    lib = C.CDLL(str(artifact), mode=os.RTLD_LOCAL | os.RTLD_NOW)
    for key in ('header_version', 'runtime_version'):
        getter = getattr(lib, 'qcb_block_'+key)
        getter.argtypes, getter.restype = [], C.c_char_p
    lib.qcb_block_abi.argtypes, lib.qcb_block_abi.restype = [], C.c_int
    lib.qcb_block_word_size.argtypes, lib.qcb_block_word_size.restype = [C.c_char_p], C.c_long
    if (lib.qcb_block_header_version(), lib.qcb_block_runtime_version(), lib.qcb_block_abi()) != (b'3.6.0', b'3.6.0', 2):
        raise RuntimeError('D13 header/runtime/bridge version mismatch')
    if any(lib.qcb_block_word_size(k.encode('ascii')) != v for k, v in _WORDS.items()):
        raise RuntimeError('D13 native word sizes mismatch')
    for name in ('solve', 'solve_rhs'):
        f = getattr(lib, 'qcb_block_'+name)
        f.argtypes = [C.c_long, C.c_long, C.c_int, C.c_void_p, C.c_size_t, C.POINTER(C.c_void_p)]
        f.restype = C.c_int
    lib.qcb_block_apply.argtypes = [C.c_long, C.c_long, C.c_long, C.c_void_p, C.c_size_t, C.POINTER(C.c_void_p)]
    lib.qcb_block_apply.restype = C.c_int
    for key, kind in (('meta', C.c_long), ('count', C.c_long), ('seconds', C.c_double)):
        getter = getattr(lib, 'qcb_block_result_'+key)
        getter.argtypes, getter.restype = [C.c_void_p, C.c_char_p], kind
    for key, kind in (('output_length', C.c_size_t), ('bytes', C.c_void_p)):
        getter = getattr(lib, 'qcb_block_result_'+key)
        getter.argtypes, getter.restype = [C.c_void_p], kind
    lib.qcb_block_result_free.argtypes, lib.qcb_block_result_free.restype = [C.c_void_p], None
    _LIBRARY = lib
    return lib


def _integer(value, minimum, maximum, label):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(label+' is outside the exact integer contract')


def _dimensions(Q, L, audit):
    _integer(Q, 1, _MAX//16, 'Q')
    _integer(L, 2, Q//2, 'L')
    if type(audit) is not bool:
        raise ValueError('audit must be a bool')
    r, B = (Q-2)//(2*L)+2, (Q+7)//8
    lengths = []
    successors = 0
    for h in (0, 1):
        T = Q-L*h
        a, b = divmod(T, 2*L)
        lengths.append(L*a+min(L, b))
        D = (T-1)//(2*L)
        successors += L*(D*T-L*D*(D+1))
    J, U = L*(2*Q-L), r*Q-r*(r-1)//2
    W = (r-1)*(Q-1)-L*(r-1)*(r-2)
    input_fields, rhs_fields = 1+(r+3)*Q, 1+r*Q+J
    output_fields = L*Q+(J+U+W+L*Q+successors if audit else 0)
    capacities = (r*Q, 3*Q, L*Q, J, U, W, successors, 2*L, 2*Q, Q+2,
        Q*(Q+1)//2, input_fields*B, rhs_fields*B, output_fields*B,
        ((Q+63)//64)*8, Q*(L+1)*8, L*Q*8, L*(L+1)*8, L*L*8,
        4096*L*Q*Q)
    if any(x < 0 or x > _MAX for x in capacities):
        raise OverflowError('D13 layout exceeds the reviewed native capacity')
    if lengths[0]+lengths[1] != Q:
        raise ArithmeticError('D13 support schedule mismatch')
    return dict(Q=Q, L=L, r=r, B=B, lengths=lengths, J=J, U=U, W=W, S=successors,
                input_fields=input_fields, rhs_fields=rhs_fields, output_fields=output_fields)


def _block_dimensions(n, L, node):
    _integer(n, 1, _MAX//16, 'n')
    _integer(L, 2, _MAX, 'L')
    _integer(node, 0, _MAX, 'node')
    k = min(L, n)
    following = n-2*L if L <= (n-1)//2 else 0
    residual = n-L if following else 0
    B = (n+7)//8
    inputs, outputs = L*n+(n-1 if k > 1 else 0), L*(k+following)
    capacities = [inputs*B, outputs*B, L*n*8, L*k*8, L*following*8,
                  ((n+63)//64)*8, 4096*L*n*n]
    if k > 1:
        capacities += [k*n*8, n*(L+1)*8, k*(L+1)*8, k*k*8]
    if following:
        capacities += [residual*L*8, L*L*8]
    if any(x < 0 or x > _MAX for x in capacities):
        raise OverflowError('D13 bare block exceeds the reviewed native capacity')
    return dict(n=n, L=L, node=node, k=k, next_bits=following,
                residual_bits=residual, B=B, input_fields=inputs, output_fields=outputs)


def _column(values, length, bits, label):
    if not isinstance(values, (list, tuple)) or len(values) != length:
        raise ValueError(label+' has an incompatible shape')
    if any(type(x) is not int or x < 0 or x.bit_length() > bits for x in values):
        raise ValueError(label+' has a noncanonical integer')


def _columns(values, count, lengths, bits, label):
    if not isinstance(values, (list, tuple)) or len(values) != count:
        raise ValueError(label+' has an incompatible column count')
    for i, column in enumerate(values):
        _column(column, lengths[i], bits[i], label)


def _read_result(lib, result, status):
    meta = {k: int(lib.qcb_block_result_meta(result, k.encode('ascii'))) for k in _META}
    counts = {k: int(lib.qcb_block_result_count(result, k.encode('ascii'))) for k in _COUNTS}
    stages = {k: float(lib.qcb_block_result_seconds(result, k.encode('ascii'))) for k in _NATIVE_STAGES}
    if meta['status'] != status or any(x < 0 for x in meta.values()) or any(x < 0 for x in counts.values()):
        raise ArithmeticError('D13 result contains unknown or invalid metadata')
    if any(not math.isfinite(x) or x < 0 for x in stages.values()):
        raise ArithmeticError('D13 result has invalid clocks')
    length = int(lib.qcb_block_result_output_length(result))
    address = lib.qcb_block_result_bytes(result)
    if status:
        if length or address:
            raise ArithmeticError('D13 failure exposed partial scientific output')
        wire = b''
    else:
        if meta['mode'] in (0, 1):
            if meta['audit'] not in (0, 1):
                raise ArithmeticError('D13 invalid result audit flag')
            layout = _dimensions(meta['Q'], meta['L'], bool(meta['audit']))
            expected = layout['output_fields']*layout['B']
            B = layout['B']
        elif meta['mode'] == 2:
            layout = _block_dimensions(meta['n'], meta['L'], meta['node'])
            if (meta['Q'], meta['audit'], meta['k'], meta['next_bits'], meta['residual_bits']) != (
                    0, 0, layout['k'], layout['next_bits'], layout['residual_bits']):
                raise ArithmeticError('D13 block metadata contradicts its schedule')
            B = layout['B']
            expected = layout['output_fields']*B
        else:
            raise ArithmeticError('D13 unknown result mode')
        if expected > _MAX or length != expected or not address or meta['field_bytes'] != B:
            raise ArithmeticError('D13 output buffer disagrees with its fixed layout')
        wire = C.string_at(address, length)
    return dict(meta=meta, counts=counts, native_stages=stages, output_length=length), wire


def _raw_request(mode, scalars, payload, *, null_data=False, null_out=False):
    """Fixed validation escape to the same C entry/result/free path.

    Scalars bypass the public mathematical checks, but must fit the ABI.
    The declared input length always equals the actual allocated bytes.
    """
    if mode not in ('solve', 'solve_rhs', 'block') or type(scalars) is not tuple or len(scalars) != 3:
        raise ValueError('invalid raw entry selector')
    if type(payload) is not bytes or type(null_data) is not bool or type(null_out) is not bool:
        raise ValueError('invalid raw byte/pointer control')
    for i, value in enumerate(scalars):
        limit = 2**31-1 if mode != 'block' and i == 2 else _MAX
        _integer(value, -limit-1, limit, 'raw scalar')
    started = perf_counter()
    t = perf_counter()
    lib = _library()
    stages = dict(load=perf_counter()-t)
    # A one-byte backing allocation keeps an empty non-NULL payload safe.
    buffer = C.create_string_buffer(payload, max(1, len(payload)))
    pointer = None if null_data else C.cast(buffer, C.c_void_p)
    result = C.c_void_p()
    status, metadata, wire, frees = None, None, b'', 0
    try:
        t = perf_counter()
        status = int(getattr(lib, 'qcb_block_'+('apply' if mode == 'block' else mode))(*scalars, pointer, len(payload),
                    None if null_out else C.byref(result)))
        stages['native_call'] = perf_counter()-t
        t = perf_counter()
        if result:
            metadata, wire = _read_result(lib, result, status)
            expected_mode = {'solve': 0, 'solve_rhs': 1, 'block': 2}[mode]
            if metadata['meta']['mode'] != expected_mode:
                raise ArithmeticError('D13 result mode disagrees with its entry')
        elif status == 0:
            raise ArithmeticError('D13 success did not return an owned result')
        stages['decode'] = perf_counter()-t
    finally:
        t = perf_counter()
        if result:
            lib.qcb_block_result_free(result)
            frees = 1
        stages['result_cleanup'] = perf_counter()-t
    elapsed = perf_counter()-started
    return dict(status=status, metadata=metadata, metadata_unavailable=metadata is None,
        output_bytes=wire, wrapper_stages=stages, result_free_calls=frees,
        total_seconds=elapsed, stage_gap_seconds=elapsed-sum(stages.values()))


def _take(wire, offset, count, B, bits):
    end = offset+count*B
    if end > len(wire):
        raise ArithmeticError('D13 truncated output segment')
    values = [int.from_bytes(wire[k:k+B], 'little') for k in range(offset, end, B)]
    if any(x.bit_length() > bits for x in values):
        raise ArithmeticError('D13 output escaped its role precision')
    return values, end


def _decode_batch(wire, layout, audit):
    Q, L, B, r = (layout[k] for k in ('Q', 'L', 'B', 'r'))
    offset, A = 0, []
    for j in range(L):
        row = []
        for _ in range(L):
            values, offset = _take(wire, offset, (Q-1-j)//L+1, B, Q)
            if any(x & ((1 << (L*z+j))-1) for z, x in enumerate(values)):
                raise ArithmeticError('D13 complete output escaped its filtration')
            row.append(values)
        A.append(row)
    arrays = None
    if audit:
        rhs, units, ratios, prefixes = [], [], [], []
        for i in range(2*L):
            bits = Q-L*(i % 2)
            values, offset = _take(wire, offset, bits, B, bits)
            rhs.append(values)
        for d in range(r):
            values, offset = _take(wire, offset, Q-d, B, Q-d)
            units.append(values)
        for d in range(r-1):
            bits = Q-2*L*d-1
            values, offset = _take(wire, offset, bits, B, bits)
            ratios.append(dict(d=d, bits=bits, coefficients=values))
        for i in range(2*L):
            values = []
            for K in range(layout['lengths'][i % 2]):
                one, offset = _take(wire, offset, 1, B, Q-L*(K//L))
                values.append(one[0])
            prefixes.append(values)
        block_next = []
        for d in range(r):
            for h in (0, 1):
                n = Q-L*h-2*L*d
                if n <= 2*L:
                    continue
                following = n-2*L
                states = []
                for _ in range(L):
                    values, offset = _take(wire, offset, following, B, following)
                    if any(x & ((1 << degree)-1) for degree, x in enumerate(values)):
                        raise ArithmeticError('D13 successor escaped filtration')
                    states.append(values)
                block_next.append(dict(d=d, h=h, n=n, k=L, next_bits=following, states=states))
        arrays = dict(normalized_rhs=rhs, units=units, unit_bits=[Q-d for d in range(r)],
                      ratios=ratios, prefixes=prefixes, block_next=block_next)
    if offset != len(wire):
        raise ArithmeticError('D13 output has trailing fields')
    return A, arrays


def _precision(layout):
    Q, L, r = (layout[k] for k in ('Q', 'L', 'r'))
    return dict(input_bits=Q, initial_bits=[Q, Q-L], prefix_lengths=layout['lengths'],
        unit_bits=[Q-d for d in range(r)], ratio_bits=[Q-2*L*d-1 for d in range(r-1)],
        blocks=[dict(d=d, h=h, n=Q-L*h-2*L*d, k=min(L, Q-L*h-2*L*d),
                     next_bits=max(Q-L*h-2*L*d-2*L, 0))
                for d in range(r) for h in (0, 1) if Q-L*h-2*L*d > 0],
        c_bits=[[Q-L*(K//L) for K in range(layout['lengths'][h])] for h in (0, 1)], output_bits=Q)


def _detail(request, started, checking, encoding, decoding, precision):
    stages = dict(check=checking, encode=encoding, **request['wrapper_stages'])
    stages['decode'] += decoding
    elapsed = perf_counter()-started
    native = request['metadata']
    return dict(status=request['status'], metadata_unavailable=native is None,
        meta=None if native is None else native['meta'],
        counts=None if native is None else native['counts'],
        native_stages=None if native is None else native['native_stages'],
        wrapper_stages=stages, precision=precision, result_free_calls=request['result_free_calls'],
        total_seconds=elapsed, stage_gap_seconds=elapsed-sum(stages.values()))


def _batch(mode, gamma, hats, columns, Q, L, audit):
    started = perf_counter()
    layout = _dimensions(Q, L, audit)
    _integer(gamma, 0, (1 << Q)-1, 'gamma')
    r, B = layout['r'], layout['B']
    _columns(hats, r, [Q]*r, [Q]*r, 'hats')
    if mode == 'solve':
        _columns(columns, 3, [Q]*3, [Q]*3, 'H')
    else:
        bits = [Q-L*(i % 2) for i in range(2*L)]
        _columns(columns, 2*L, bits, bits, 'normalized RHS')
    checked = perf_counter()
    payload = b''.join(x.to_bytes(B, 'little') for group in ([gamma], *hats, *columns) for x in group)
    encoded = perf_counter()
    request = _raw_request(mode, (Q, L, int(audit)), payload)
    t = perf_counter()
    arrays, A = None, None
    if request['status'] == 0:
        m = request['metadata']['meta']
        if (m['Q'], m['L'], m['audit'], m['n'], m['node'], m['k'], m['next_bits'], m['residual_bits']) != (Q, L, int(audit), 0, 0, 0, 0, 0):
            raise ArithmeticError('D13 native result task mismatch')
        A, arrays = _decode_batch(request['output_bytes'], layout, audit)
    decoding = perf_counter()-t
    detail = _detail(request, started, checked-started, encoded-checked, decoding, _precision(layout))
    if request['status']:
        raise NativeBlockError(request['status'], detail)
    if audit:
        detail['audit_arrays'] = arrays
    # Metadata assembly is part of the enclosing bridge cost, including gap.
    detail['total_seconds'] = perf_counter()-started
    detail['stage_gap_seconds'] = detail['total_seconds']-sum(detail['wrapper_stages'].values())
    return A, detail


def second_transform(gamma, hats, H, Q, L, *, audit=False):
    return _batch('solve', gamma, hats, H, Q, L, audit)


def solve_rhs(gamma, hats, normalized_rhs, Q, L, *, audit=False):
    return _batch('solve_rhs', gamma, hats, normalized_rhs, Q, L, audit)


def block(values, node, rho, n, L):
    """Apply the same full block core used by the production constructor."""
    started = perf_counter()
    layout = _block_dimensions(n, L, node)
    k, following, B = (layout[name] for name in ('k', 'next_bits', 'B'))
    _columns(values, L, [n]*L, [n]*L, 'block states')
    if k == 1:
        if rho is not None:
            raise ValueError('one-coefficient block must not receive rho')
        rho_values = []
    else:
        _column(rho, n-1, n-1, 'block rho')
        rho_values = rho
    checked = perf_counter()
    payload = b''.join(x.to_bytes(B, 'little') for group in (*values, rho_values) for x in group)
    encoded = perf_counter()
    request = _raw_request('block', (n, L, node), payload)
    t = perf_counter()
    a, successor = None, None
    if request['status'] == 0:
        m = request['metadata']['meta']
        if (m['n'], m['L'], m['node'], m['k'], m['next_bits'], m['residual_bits'],
                m['Q'], m['audit']) != (n, L, node, k, following, layout['residual_bits'], 0, 0):
            raise ArithmeticError('D13 bare block result task mismatch')
        offset, a = 0, []
        for _ in range(L):
            row = []
            for j in range(k):
                item, offset = _take(request['output_bytes'], offset, 1, B, n-j)
                row.append(item[0])
            a.append(row)
        if following:
            successor = []
            for _ in range(L):
                row, offset = _take(request['output_bytes'], offset, following, B, following)
                if any(x & ((1 << degree)-1) for degree, x in enumerate(row)):
                    raise ArithmeticError('D13 bare successor escaped filtration')
                successor.append(row)
        if offset != len(request['output_bytes']):
            raise ArithmeticError('D13 block trailing output fields')
    decoding = perf_counter()-t
    detail = _detail(request, started, checked-started, encoded-checked, decoding,
                     dict(input_bits=n, rho_bits=n-1 if k > 1 else None,
                          coefficient_bits=[n-j for j in range(k)],
                          next_bits=following or None, node=node))
    if request['status']:
        raise NativeBlockError(request['status'], detail)
    detail['total_seconds'] = perf_counter()-started
    detail['stage_gap_seconds'] = detail['total_seconds']-sum(detail['wrapper_stages'].values())
    return a, successor, detail

