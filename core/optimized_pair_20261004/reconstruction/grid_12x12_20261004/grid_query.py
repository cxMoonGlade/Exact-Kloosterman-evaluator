"""Grid-capacity-extended D15 parameter/specialization and resident ordered readout.

Callers certify the actual original f/a with common.validate_input. Both legs
use the same support envelope, frozen D14 new parameter algorithm, and P/Q
precision. Native owners never cross query boundaries. Import does no math.
"""
from time import perf_counter

import grid_parameter_query as _old
import grid_interop as bridge

IMPLEMENTATION = 'resident_ordered_readout_grid12_v1'
READOUTS = ('old', 'new')


def _precheck(acquired, method, N, L, readout, keep_debug):
    if type(readout) is not str or readout not in READOUTS:
        raise ValueError('unknown registered readout')
    if type(keep_debug) is not bool:
        raise ValueError('keep_debug must be a bool')
    if method not in ('ah', 'ordinary') or not isinstance(acquired, dict):
        raise ValueError('invalid acquisition or method')
    P = acquired.get('P', acquired.get('output_P'))
    # d is not precomputed or guessed: 1 is only a harmless scalar-envelope
    # witness here. Both actual algorithms subsequently obtain and check d.
    bridge.checked_input(N, L, 1, P, keep_debug)
    return P


def _matrix_from_export(R, exported):
    meta = exported['meta']
    if meta['parts'] != 2:
        raise ArithmeticError('D15 matrix export lost its Gaussian component')
    return [[(R.ctx(pair[0]), R.ctx(pair[1])) for pair in row]
            for row in exported['coefficients']]


def query(acquired, method, N, L, f, a, readout='old', keep_debug=False):
    """Compute the exact K_L(a), paying all conversion/load/cleanup locally."""
    started = perf_counter()
    P = _precheck(acquired, method, N, L, readout, keep_debug)
    prechecked = perf_counter()
    if readout == 'old':
        row = _old.query(acquired, method, N, L, f, a,
                         parameter_method='new', keep_debug=keep_debug)
        row['readout'] = 'old'
        overhead = prechecked-started
        row['parameter_detail']['input_check'] += overhead
        row['stages']['parameter'] += overhead
        # Retain the unchanged old operation statistics. The small wrapper tail
        # is charged as readout work, not silently removed from the local total.
        elapsed = perf_counter()-started
        row['stages']['binary_product'] = elapsed-row['stages']['parameter']-row['stages']['specialization']
        row['stages']['total'] = elapsed
        return row

    # These checks/preparation follow the frozen D14 query body. The only
    # changed mathematical code begins after specialization.
    _old.parameter.checked_input(N, f, a)
    if acquired.get('N') != N or acquired.get('L') != L:
        raise ValueError('acquisition metadata must use the actual N and L')
    Q = acquired.get('Q') if method == 'ah' else P
    if type(Q) is not int or Q != (P+L-1 if method == 'ah' else P):
        raise ValueError('AH row division requires the registered guard')
    checked = perf_counter()
    fa, d, algorithm = _old.parameter.parameter_polynomial(N, f, a)
    if (type(d) is not int or d < 1 or d > N or N % d or
            not isinstance(fa, list) or len(fa) != d+1 or
            any(type(c) is not int or c not in (0, 1) for c in fa) or fa[-1] != 1):
        raise ArithmeticError('parameter construction returned an invalid polynomial or degree')
    polynomial_done = perf_counter()
    parent = _old._tq.ForwardPresentation(fa, Q)
    prepared = perf_counter()
    B = _old._tq.common.specialize(acquired, method, parent, L, P)
    specialized = perf_counter()
    R = parent.ring(P)
    gamma = acquired['gamma']
    if type(gamma) is not int:
        raise ValueError('gamma must be an exact Python integer')
    bridge.checked_input(N, L, d, P, keep_debug)
    A = [[R.scale(value, pow(gamma, i, R.mod)) for value in row]
         for i, row in enumerate(B)]
    # Public getters create independent, canonical P-bit wire inputs. No
    # Python-flint object layout or C pointer is borrowed.
    g = [int(R.g[k]) for k in range(d+1)]
    coefficients = [[[[int(poly[k]) for k in range(d)] for poly in pair]
                     for pair in row] for row in A]
    owner = bridge.ResidentReadout(N, L, d, P, g, coefficients, audit=keep_debug)
    trace = product = pair = debug_steps = None
    try:
        owner.run()
        exported = owner.export(bridge.TRACE)
        trace_matrix = _matrix_from_export(R, exported)
        if len(trace_matrix) != 1 or len(trace_matrix[0]) != 1:
            raise ArithmeticError('D15 trace did not use its exact 1 by 1 output shape')
        trace = trace_matrix[0][0]
        if keep_debug:
            before_close = owner.statistics_snapshot()
            debug_steps = dict(initial_A=owner.export(bridge.INITIAL_A), steps=[])
            for index, meta in enumerate(before_close['native']['steps']):
                entry = dict(meta=meta, left=owner.export(bridge.LEFT, index),
                             right=owner.export(bridge.RIGHT, index),
                             product=owner.export(bridge.PRODUCT, index),
                             jump_before=None, jump_after=None)
                if meta['phase'] == 0:
                    entry['jump_before'] = owner.export(bridge.JUMP_BEFORE, index)
                    entry['jump_after'] = owner.export(bridge.JUMP_AFTER, index)
                debug_steps['steps'].append(entry)
            if debug_steps['steps']:
                last = debug_steps['steps'][-1]
                pair = (_matrix_from_export(R, last['left']), _matrix_from_export(R, last['right']))
            else:
                product = _matrix_from_export(R, owner.export(bridge.FINAL_MATRIX))
        value = _old._tq.recover_trace(R, trace, L)
    finally:
        owner.close()
    detail = owner.statistics_snapshot()
    statistics = dict(degree=d, bits=P, rank=L, subfield_power=N//d,
                      parent_xi_constructed=parent.xi is not None,
                      native=detail['native'])
    row = dict(implementation=IMPLEMENTATION, mode='combined', readout='new', a=a, value=value,
               degree=d, fa=fa, P=P, Q=Q, operations=statistics, parameter_method='new',
               parameter_detail=dict(input_check=checked-started,
                   polynomial=polynomial_done-checked, presentation=prepared-polynomial_done,
                   algorithm=algorithm), resident_detail=detail)
    if keep_debug:
        row['_debug'] = dict(parent=parent, B=B, trace=trace, final_matrix=product,
                             terminal_U=pair[0] if pair is not None else None,
                             terminal_V=pair[1] if pair is not None else None,
                             resident_steps=debug_steps)
    # Destruction of native objects happened in close; dropping temporary host
    # arrays and polynomial references is also within binary_product.
    parent = B = R = A = g = coefficients = owner = trace = product = pair = debug_steps = None
    exported = trace_matrix = None
    ended = perf_counter()
    row['stages'] = dict(parameter=prepared-started, specialization=specialized-prepared,
                         binary_product=ended-specialized, total=ended-started)
    return row
