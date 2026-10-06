"""D13 complete cold AH acquisition, at the shared sufficient precision.

The first guard and scalar constructor are unchanged. The second
operation uses the shared block core in an owned binary C bridge. No mathematical cache is kept.
"""
import importlib
from pathlib import Path
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
D10 = HERE.parent/'carrier_coverage'
for _path in (HERE, D10):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
_fair = importlib.import_module('fair_acquire')
_bridge = importlib.import_module('block_interop')
if (Path(_fair.__file__).resolve() != D10/'fair_acquire.py'
        or Path(_bridge.__file__).resolve() != HERE/'block_interop.py'):
    raise RuntimeError('unexpected D13 precision/bridge module source')


def _modules():
    ah = importlib.import_module('ah_native')
    guard = importlib.import_module('guard_interop')
    if (Path(ah.__file__).resolve() != ROOT/'code/ah_native.py'
            or Path(guard.__file__).resolve() != HERE.parent/'precision/guard_interop.py'):
        raise RuntimeError('unexpected D13 scalar/first-guard module source')
    import flint
    if (flint.__version__, flint.__FLINT_VERSION__, flint.ctx.threads) != ('0.9.0', '3.6.0', 1):
        raise RuntimeError('D13 requires the reviewed single-thread FLINT version')
    return ah, guard


def acquire(N, L, P=None, *, audit=False):
    started = perf_counter()
    plan = _fair.precision_plan(N, L, P)
    actual_P, Q = plan['requested_P'], plan['requested_P']+L-1
    layout = _bridge._dimensions(Q, L, audit)
    r = layout['r']
    if Q*(r+3) > sys.maxsize or 2*(Q-1)*(r+3) > sys.maxsize:
        raise OverflowError('D13 first-guard layout exceeds the native ABI')
    checked = perf_counter()
    ah, guard = _modules()
    loaded = perf_counter()
    stages = dict(input_checks=checked-started, backend_loading=loaded-checked)
    mask = (1 << Q)-1

    t = perf_counter()
    gamma, z, phi_values, b = ah._scalar_data(Q)
    if len(phi_values) != Q or len(b) < 2*Q:
        raise ValueError('incompatible cold scalar-constructor output')
    stages['scalar_construction'] = perf_counter()-t

    t = perf_counter()
    gamma_powers = [1]
    for _ in range(Q):
        gamma_powers.append((gamma_powers[-1]*gamma) & mask)
    H = [[0 if 2*k+offset < 0 else
          (gamma_powers[k+(offset == 1)]*b[2*k+offset]) & mask
          for k in range(Q)] for offset in (-1, 0, 1)]
    del b, z, gamma_powers
    stages['first_rhs_construction'] = perf_counter()-t

    t = perf_counter()
    first_solved, first_detail = guard.first_transform(phi_values, H, Q, r)
    _bridge._columns(first_solved, r+3, [Q]*(r+3), [Q]*(r+3), 'first guard output')
    hats, h_values = first_solved[:r], first_solved[r:]
    if not audit:
        del first_solved
    del phi_values, H
    stages['first_shared_transform'] = perf_counter()-t

    t = perf_counter()
    coefficients, second_detail = _bridge.second_transform(gamma, hats, h_values, Q, L, audit=audit)
    del hats, h_values
    stages['second_block'] = perf_counter()-t
    computed = perf_counter()

    # This constructor directly uses Q,L. It never calls acquire(N=1).
    plan.update(backend='ah_block_second', backend_module=__name__,
                backend_source=str(Path(__file__).resolve()), backend_call_N=None,
                backend_dimension='Q,L', coefficient_constructor_Q=Q)
    result = dict(N=N, L=L, P=actual_P, Q=Q, gamma=gamma, coefficients=coefficients,
        stages=stages,
        counts=dict(first_solve_rhs=r+3, first_input_matrix_slots=0,
            first_rhs_slots=3*Q, first_output_slots=Q*(r+3),
            native_factorizations=0, second_solve_rhs=2*L, second_input_matrix_slots=0,
            second_rhs_slots=L*(2*Q-L), partial_c_items=L*Q, supported_output_slots=L*Q),
        first_transform_detail=first_detail, second_block_detail=second_detail,
        precision=dict(output=actual_P, acquisition=Q,
            intermediate_solution='visible prefixes only, per-item Q-Ld bits'),
        diagnostics=dict(enabled=audit, full_system_residuals=0),
        native=dict(python_flint='0.9.0', flint='3.6.0', second_bridge_abi=2),
        layout='coefficients[target_row][source_column][degree]',
        algorithm='guard_first_transform_then_shared_block_elimination',
        recovery_precision=plan)
    if audit:
        result['audit_arrays'] = dict(second_detail['audit_arrays'], first_solved=first_solved)
    adapter = dict(input_checks=checked-started, backend_loading=loaded-checked,
                   backend_call=computed-loaded)
    plan['adapter_stages'] = adapter
    finished = perf_counter()
    adapter.update(metadata=finished-computed, total=finished-started)
    result.update(acquisition_seconds=finished-started,
                  stage_gap_seconds=finished-started-sum(stages.values()))
    return result

