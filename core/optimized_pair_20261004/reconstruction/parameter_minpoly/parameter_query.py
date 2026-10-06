"""D14 shared parameter construction with the frozen D11 combined readout.

This internal function requires an already certified irreducible binary f and
nonzero original-basis a. Registered callers pay common.validate_input before
entry. Shape checks here do not certify irreducibility. No scientific operation
or bridge loading occurs on import; there is no parameter or query cache.
"""
from time import perf_counter

import parameter
import terminal_query as _tq


IMPLEMENTATION = "shared_parameter_readout_d14_v1"
PARAMETER_METHODS = ("old", "new")


def query(acquired, method, N, L, f, a, parameter_method="old", keep_debug=False):
    """Compute the same exact integer from a paid, certified original input.

Both parameter methods use this wrapper. Only acquisition of fa and its actual
degree changes; the high-precision presentation, specialization, ordered trace
and integer recovery are the frozen D11 operations. Debug objects are retained
only when explicitly requested by the registered validator.
    """
    started = perf_counter()
    if type(parameter_method) is not str or parameter_method not in PARAMETER_METHODS:
        raise ValueError("unknown registered parameter method")
    if type(keep_debug) is not bool:
        raise ValueError("keep_debug must be a boolean")
    if method not in ("ah", "ordinary") or not isinstance(acquired, dict):
        raise ValueError("invalid acquisition or method")
    parameter.checked_input(N, f, a)
    if type(L) is not int or L < 2:
        raise ValueError("expected L>=2")
    if acquired.get("N") != N or acquired.get("L") != L:
        raise ValueError("acquisition metadata must use the actual N and L")
    P = acquired.get("P", acquired.get("output_P"))
    if type(P) is not int or P < 1:
        raise ValueError("acquisition must specify its actual output precision")
    Q = acquired.get("Q") if method == "ah" else P
    if type(Q) is not int or Q != (P + L - 1 if method == "ah" else P):
        raise ValueError("AH row division requires the registered guard")
    checked = perf_counter()

    if parameter_method == "old":
        fa, orbit = _tq.common.parameter_polynomial(a, f)
        d = len(orbit)
        # Dropping the actual orbit is charged to the polynomial stage.
        orbit = None
        algorithm = dict(method="frozen_conjugate_product_tree", N=N, degree=d,
                         counts=None, seconds=None)
    else:
        fa, d, algorithm = parameter.parameter_polynomial(N, f, a)
    if (type(d) is not int or d < 1 or d > N or N % d or
            not isinstance(fa, list) or len(fa) != d + 1 or
            any(type(c) is not int or c not in (0, 1) for c in fa) or fa[-1] != 1):
        raise ArithmeticError("parameter construction returned an invalid polynomial or degree")
    polynomial_done = perf_counter()

    parent = _tq.ForwardPresentation(fa, Q)
    prepared = perf_counter()
    B = _tq.common.specialize(acquired, method, parent, L, P)
    specialized = perf_counter()
    trace, product, statistics, pair = _tq.ordered_trace(
        parent, B, acquired["gamma"], N, L, P, "combined", keep_debug)
    value = _tq.recover_trace(parent.ring(P), trace, L)
    row = dict(implementation=IMPLEMENTATION, mode="combined", a=a, value=value,
               degree=d, fa=fa, P=P, Q=Q, operations=statistics,
               parameter_method=parameter_method,
               parameter_detail=dict(input_check=checked-started,
                   polynomial=polynomial_done-checked,
                   presentation=prepared-polynomial_done, algorithm=algorithm))
    if keep_debug:
        row["_debug"] = dict(
            parent=parent, B=B, trace=trace, final_matrix=product,
            terminal_U=pair[0] if pair is not None else None,
            terminal_V=pair[1] if pair is not None else None,
        )
    parent = B = trace = product = pair = None
    ended = perf_counter()
    row["stages"] = dict(parameter=prepared-started,
                         specialization=specialized-prepared,
                         binary_product=ended-specialized, total=ended-started)
    return row
