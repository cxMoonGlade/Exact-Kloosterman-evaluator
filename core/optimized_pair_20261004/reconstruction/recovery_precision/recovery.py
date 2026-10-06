"""D9: explicit shared recovery precision around four frozen acquisitions.

Only module loading is cached. Every acquisition recomputes its coefficients;
the internal call N=1 selects a precision-independent coefficient primitive.
Queries must still receive the actual N and original field/parameter input.
Import performs no mathematical computation or native-library loading.
"""
import importlib
import importlib.util
from pathlib import Path
import sys
from time import perf_counter


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
BACKENDS = ('ah_native', 'ah_guard', 'ah_output', 'ordinary')
POLICY = 'finite_fourier_shared_sufficient_v1'
_GUARD_NAME = '_qcb_d9_guard_acquisition'

# These are import search locations, not imported modules or loaded libraries.
for _path in (ROOT/'code', HERE.parent/'precision',
              HERE.parent/'output_precision', HERE.parent/'first_transform'):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))


def precision_plan(N, L, P=None):
    """Validate actual task integers and return the sufficient precision policy."""
    if type(N) is not int or type(L) is not int or N < 1 or L < 2:
        raise ValueError('positive integer N and integer L>=2 required')
    old = N*(L-1)+2
    shared = max(L+1, (N*L+1)//2+1)
    effective = shared if P is None else P
    if type(effective) is not int or effective < shared:
        raise ValueError('P is below the shared sufficient recovery precision')
    return dict(N=N, L=L, target_N=N, backend_call_N=1,
                P0=old, Pshared=shared, requested_P=effective,
                backend_minimum_output_P=L+1, policy=POLICY,
                selection='default_shared' if P is None else 'explicit')


def _check_backend(backend):
    if type(backend) is not str or backend not in BACKENDS:
        raise ValueError('unknown acquisition backend')


def _checked_module(name, path):
    module = importlib.import_module(name)
    if Path(module.__file__).resolve() != path.resolve():
        raise RuntimeError('unexpected module source for '+name)
    return module


def _guard_module():
    """Load an owned module instance without rebinding the frozen baseline."""
    guard = _checked_module('guard_interop', HERE.parent/'precision/guard_interop.py')
    path = HERE.parent/'first_transform/first_acquisition.py'
    if _GUARD_NAME in sys.modules:
        result = sys.modules[_GUARD_NAME]
        if (Path(result.__file__).resolve() != path.resolve() or
                result.first_transform is not guard.first_transform):
            raise RuntimeError('changed owned guard acquisition binding')
        return result
    sentinel = object()
    previous = sys.modules.get('first_interop', sentinel)
    baseline = sys.modules.get('first_acquisition')
    previous_function = None if baseline is None else baseline.first_transform
    spec = importlib.util.spec_from_file_location(_GUARD_NAME, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[_GUARD_NAME] = result
    try:
        sys.modules['first_interop'] = guard
        spec.loader.exec_module(result)
        if result.first_transform is not guard.first_transform:
            raise RuntimeError('owned acquisition did not capture the guard')
    except BaseException:
        sys.modules.pop(_GUARD_NAME, None)
        raise
    finally:
        if previous is sentinel:
            sys.modules.pop('first_interop', None)
        else:
            sys.modules['first_interop'] = previous
        if sys.modules.get('first_interop', sentinel) is not previous:
            raise RuntimeError('first_interop module restoration failed')
        if baseline is not None and baseline.first_transform is not previous_function:
            raise RuntimeError('baseline acquisition binding changed')
    return result


def backend_module(backend):
    """Load only the selected frozen coefficient primitive, without executing it."""
    _check_backend(backend)
    if backend == 'ah_guard':
        return _guard_module()
    name, path = {
        'ah_native': ('ah_native', ROOT/'code/ah_native.py'),
        'ah_output': ('output_acquisition', HERE.parent/'output_precision/output_acquisition.py'),
        'ordinary': ('ordinary', ROOT/'code/ordinary.py'),
    }[backend]
    result = _checked_module(name, path)
    if backend == 'ah_output':
        guard = _checked_module('guard_interop', HERE.parent/'precision/guard_interop.py')
        if result.first_transform is not guard.first_transform:
            raise RuntimeError('changed output-acquisition guard binding')
    return result


def _early_guard_dimensions(plan):
    """Mirror the frozen Python wrapper's three long bounds before scalars.

    The C primitive retains its additional tree/allocation-size checks. This
    early check is not an allocation-feasibility test or a complete C guard.
    """
    import ctypes
    Q, L = plan['requested_P']+plan['L']-1, plan['L']
    r = (Q+L-2)//L+1
    maximum = (1 << (8*ctypes.sizeof(ctypes.c_long)-1))-1
    if Q > maximum//16 or Q*(r+3) > maximum or 2*(Q-1)*(r+3) > maximum:
        raise ValueError('D5 dimensions exceed the native ABI')


def acquire(backend, N, L, P=None, *, audit=False):
    """Recompute coefficients at explicit P, defaulting to shared precision.

    ordinary has no audit operation and rejects audit=True. AH retains its
    frozen audit operation. All adapter work is included in adapter_stages;
    the frozen backend's own stage/elapsed fields are left unchanged.
    """
    started = perf_counter()
    _check_backend(backend)
    plan = precision_plan(N, L, P)
    if type(audit) is not bool:
        raise ValueError('audit must be a bool')
    if backend == 'ordinary' and audit:
        raise ValueError('ordinary acquisition does not define audit=True')
    if backend == 'ah_guard':
        _early_guard_dimensions(plan)
    checked = perf_counter()
    module = backend_module(backend)
    loaded = perf_counter()
    actual_P = plan['requested_P']
    if backend == 'ordinary':
        result = module.acquire(1, L, P=actual_P)
    else:
        result = module.acquire(1, L, P=actual_P, audit=audit)
    called = perf_counter()
    if backend == 'ordinary':
        if result['output_P'] != actual_P or result['precision']['P'] != actual_P:
            raise ArithmeticError('ordinary backend returned the wrong precision')
        result['precision']['N'] = N
        result['precision']['minimum_output_P'] = plan['Pshared']
    elif result['P'] != actual_P or result['Q'] != actual_P+L-1:
        raise ArithmeticError('AH backend returned the wrong precision')
    result['N'], result['L'] = N, L
    plan.update(backend=backend, backend_module=module.__name__,
                backend_source=str(Path(module.__file__).resolve()))
    result['recovery_precision'] = plan
    stages = dict(input_checks=checked-started, backend_loading=loaded-checked,
                  backend_call=called-loaded)
    plan['adapter_stages'] = stages
    finished = perf_counter()
    stages.update(metadata=finished-called, total=finished-started)
    return result
