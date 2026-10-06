"""D10: one acquisition adapter for four frozen backends and ordinary_int.

The D9 coefficient-independent precision plan and selective module loader
are reused. Every backend is acquired directly through the same timed body;
no backend goes through another acquire wrapper or reuses computed arrays.
Import loads only standard-library code and frozen D9 adapter definitions.
"""
import importlib
from pathlib import Path
import sys
from time import perf_counter


HERE = Path(__file__).resolve().parent
D9 = HERE.parent / "recovery_precision"
for _path in (HERE, D9):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
_legacy = importlib.import_module("recovery")
if Path(_legacy.__file__).resolve() != D9 / "recovery.py":
    raise RuntimeError("unexpected D9 adapter source")

BACKENDS = ("ah_native", "ah_guard", "ah_output", "ordinary", "ordinary_int")
ORDINARY = ("ordinary", "ordinary_int")


def precision_plan(N, L, P=None):
    """Use the reviewed sufficient precision for the actual task N and L."""
    return _legacy.precision_plan(N, L, P)


def _check_backend(backend):
    if type(backend) is not str or backend not in BACKENDS:
        raise ValueError("unknown acquisition backend")


def backend_module(backend):
    """Load only the chosen module, with its actual source path checked."""
    _check_backend(backend)
    if backend == "ordinary_int":
        return _legacy._checked_module("ordinary_int", HERE / "ordinary_int.py")
    return _legacy.backend_module(backend)


def acquire(backend, N, L, P=None, *, audit=False):
    """Recompute complete coefficients and restore true-N task metadata.

Both ordinary implementations keep the same initializer and full precision
schedule; neither defines audit=True. The actual backend's seconds/stages,
counts/counters and implementation identity are preserved without relabeling.
Queries must use the original N/f/a, never the internal coefficient call N=1.
"""
    started = perf_counter()
    _check_backend(backend)
    plan = precision_plan(N, L, P)
    if type(audit) is not bool:
        raise ValueError("audit must be a bool")
    ordinary = backend in ORDINARY
    if ordinary and audit:
        raise ValueError("ordinary acquisition does not define audit=True")
    if backend == "ah_guard":
        _legacy._early_guard_dimensions(plan)
    checked = perf_counter()
    module = backend_module(backend)
    loaded = perf_counter()
    actual_P = plan["requested_P"]
    if ordinary:
        result = module.acquire(1, L, P=actual_P)
    else:
        result = module.acquire(1, L, P=actual_P, audit=audit)
    called = perf_counter()
    if ordinary:
        if result["output_P"] != actual_P or result["precision"]["P"] != actual_P:
            raise ArithmeticError("ordinary backend returned the wrong precision")
        result["precision"]["N"] = N
        result["precision"]["minimum_output_P"] = plan["Pshared"]
    elif result["P"] != actual_P or result["Q"] != actual_P + L - 1:
        raise ArithmeticError("AH backend returned the wrong precision")
    result["N"], result["L"] = N, L
    plan.update(backend=backend, backend_module=module.__name__,
                backend_source=str(Path(module.__file__).resolve()))
    result["recovery_precision"] = plan
    stages = dict(input_checks=checked-started, backend_loading=loaded-checked,
                  backend_call=called-loaded)
    plan["adapter_stages"] = stages
    finished = perf_counter()
    stages.update(metadata=finished-called, total=finished-started)
    return result
