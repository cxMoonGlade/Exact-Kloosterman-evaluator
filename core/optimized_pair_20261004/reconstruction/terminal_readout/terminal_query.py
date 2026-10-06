"""D11: three exact, shared changes to the frozen D6 binary readout.

No scientific operation or native bridge load occurs on import. The five
modes use the same observation and accounting skeleton. The packed kernel,
parameter arithmetic, specialization, and later modular compositions remain
the frozen implementations. No reference or diagnostic module is imported.
"""
from struct import calcsize
import sys
from time import perf_counter

from flint import fmpz_mod_poly

import common
from packed import PackedMultiplier


IMPLEMENTATION = "terminal_shared_readout_d11_v1"
MODES = ("base", "lazy_tau", "first_inflate", "terminal_trace", "combined")
_LONG_MAX = (1 << 63) - 1
_SIZE_MAX = (1 << 64) - 1
_COUNTERS = (
    "calls", "successful_calls", "failed_calls", "zero_calls", "real_calls",
    "one_imag_calls", "gaussian_calls", "native_products", "mod_g_calls",
    "input_coefficients", "packed_entries", "packed_coefficients",
    "decoded_coefficients",
)


def _mode(mode):
    if type(mode) is not str or mode not in MODES:
        raise ValueError("unknown D11 readout mode")
    return mode


def checked_inflate_length(length):
    """Check the actual public inflate(2) output length without allocating it.

    The registered native platform is LP64. This bounds explicit coefficient
    storage, not FLINT's entire temporary workspace; process limits still apply.
    Zero and constant polynomials take the public C implementation's copy path.
    """
    if (calcsize("l"), calcsize("L"), calcsize("P")) != (8, 8, 8):
        raise RuntimeError("D11 inflate contract requires the registered LP64 ABI")
    if type(length) is not int or length < 0:
        raise ValueError("polynomial length must be a nonnegative Python integer")
    target = length if length <= 1 else 2 * (length - 1) + 1
    if target > _LONG_MAX or target > min(_SIZE_MAX, sys.maxsize) // 8:
        raise OverflowError("inflated coefficient length exceeds the size contract")
    return target


def _check_component(R, polynomial):
    if not isinstance(polynomial, fmpz_mod_poly):
        raise ValueError("expected a public modular polynomial")
    if int(polynomial.modulus()) != R.mod:
        raise ValueError("polynomial coefficient context mismatch")
    if polynomial.degree() >= R.d:
        raise ValueError("expected a reduced polynomial of degree below d")


def inflate_once(R, polynomial):
    """Return p(X**2) mod g, retaining each coefficient and fixing Gaussian i."""
    _check_component(R, polynomial)
    checked_inflate_length(len(polynomial))
    return polynomial.inflate(2) % R.g


class ForwardPresentation(common.Presentation):
    """The same paid forward presentation, with inherited lazy true-tau use."""

    def __init__(self, fa, e):
        self.g = common.teich_modulus(fa, e)
        self.e, self.cache, self.xi = e, {}, None


def _square_shape(matrix):
    if not isinstance(matrix, (list, tuple)) or not matrix:
        raise ValueError("expected a nonempty square matrix")
    n = len(matrix)
    if any(not isinstance(row, (list, tuple)) or len(row) != n for row in matrix):
        raise ValueError("expected a square matrix")
    if n > min(_LONG_MAX, sys.maxsize) // n:
        raise OverflowError("matrix entry count exceeds the size contract")
    if n * n > min(_LONG_MAX, sys.maxsize) // calcsize("P"):
        raise OverflowError("matrix entry storage exceeds the size contract")
    return n


def terminal_trace(multiplier, U, V):
    """Compute the whole Gaussian trace via one 1 x L**2 by L**2 x 1 product.

    The frozen multiplier validates all component contexts, canonical residues,
    degrees, and the new inner dimension before packing. There is no conjugation.
    """
    L = _square_shape(U)
    if _square_shape(V) != L:
        raise ValueError("terminal factors have incompatible dimensions")
    left = [[U[i][j] for i in range(L) for j in range(L)]]
    right = [[V[j][i]] for i in range(L) for j in range(L)]
    return multiplier.multiply(left, right)[0][0]


def recover_trace(R, trace, L):
    """Check the complete element before the unchanged signed integer recovery."""
    if type(L) is not int or L < 2:
        raise ValueError("rank must be a Python integer at least two")
    if not isinstance(trace, (list, tuple)) or len(trace) != 2:
        raise ValueError("expected a Gaussian polynomial pair")
    for component in trace:
        _check_component(R, component)
    if trace[1] or trace[0].degree() > 0:
        raise ArithmeticError("complete trace is not a real base-ring scalar")
    return common.centered((-1) ** (L - 1) * int(trace[0][0]), R.e)


def _product_count(exponent):
    return exponent.bit_length() + exponent.bit_count() - 2


def ordered_trace(parent, B, gamma, N, L, P, mode="base", keep_debug=False):
    """Return (whole_trace, final_matrix_or_None, statistics, final_pair_or_None).

    The final pair is retained only when keep_debug=True and a multiplication
    actually exists. Terminal N>1 modes never construct the final matrix.
    N=1 already has A, so returning that existing matrix does not fabricate one.
    """
    _mode(mode)
    if type(keep_debug) is not bool:
        raise ValueError("keep_debug must be a boolean")
    if type(N) is not int or N < 1 or type(L) is not int or L < 2:
        raise ValueError("expected N>=1 and L>=2")
    if type(P) is not int or P < 1 or type(gamma) is not int:
        raise ValueError("invalid precision or scalar")
    R = parent.ring(P)
    d = R.d
    if N % d or _square_shape(B) != L:
        raise ValueError("incompatible degree or matrix dimensions")
    e = N // d
    use_inflate = mode in ("first_inflate", "combined")
    use_terminal = mode in ("terminal_trace", "combined")
    expected_products = _product_count(d) + _product_count(e)
    multiplier = PackedMultiplier(R)
    statistics = dict(
        degree=d, bits=P, rank=L, subfield_power=e,
        frobenius_matrix_products=0, subfield_matrix_products=0,
        full_matrix_products=0, terminal_trace_products=0,
        component_compositions=0, argument_compositions=0,
        argument_squares=0, doublings=0, appends=0, batch_calls=0,
        inflation_components=0, inflation_jump_calls=0,
        composer_created=0, composer_closed=True, packed_calls=[],
        parent_xi_constructed=parent.xi is not None,
    )
    A = [[R.scale(value, pow(gamma, i, R.mod)) for value in row]
         for i, row in enumerate(B)]
    product = A
    composer = None
    final_trace = None
    final_pair = None
    products_done = 0

    def multiply(U, V):
        nonlocal products_done, final_trace, final_pair
        products_done += 1
        last = products_done == expected_products
        if products_done > expected_products:
            raise ArithmeticError("binary product schedule exceeded its fixed length")
        if last and keep_debug:
            final_pair = U, V
        terminal = use_terminal and last
        kind = "terminal" if terminal else "full"
        shape = [1, L * L, 1] if terminal else [L, L, L]
        has_ai = any(bool(value[1]) for row in U for value in row)
        has_bi = any(bool(value[1]) for row in V for value in row)
        before = multiplier.statistics_snapshot()
        try:
            if terminal:
                final_trace = terminal_trace(multiplier, U, V)
                statistics["terminal_trace_products"] += 1
                return None
            result = multiplier.multiply(U, V)
            statistics["full_matrix_products"] += 1
            return result
        finally:
            after = multiplier.statistics_snapshot()
            guard = (shape[1] * d - 1).bit_length() + 1
            field_bits = 8 * ((2 * P + guard + 7) // 8)
            statistics["packed_calls"].append(dict(
                kind=kind, shape=shape, field_bits=field_bits,
                has_left_imag=has_ai, has_right_imag=has_bi,
                deltas={name: after[name] - before[name] for name in _COUNTERS},
            ))

    def compose(matrix, argument, include_jump, first):
        nonlocal composer
        if use_inflate and first:
            values = []
            for row in matrix:
                new_row = []
                for value in row:
                    parts = []
                    for polynomial in value:
                        if polynomial:
                            parts.append(inflate_once(R, polynomial))
                            statistics["inflation_components"] += 1
                        else:
                            parts.append(polynomial)
                    new_row.append(tuple(parts))
                values.append(new_row)
            if include_jump:
                argument = inflate_once(R, argument)
                statistics["inflation_jump_calls"] += 1
            return values, argument

        if composer is None:
            # Both module import and owned native context acquisition are paid
            # here, only if a genuine general composition is required.
            from interop_flint import BatchComposer
            composer = BatchComposer(R)
            statistics["composer_created"] = 1
            statistics["composer_closed"] = False
        positions, outer = [], []
        for i, row in enumerate(matrix):
            for j, value in enumerate(row):
                for part in range(2):
                    if value[part]:
                        positions.append((i, j, part))
                        outer.append(value[part])
        statistics["component_compositions"] += len(outer)
        if include_jump:
            outer.append(argument)
            statistics["argument_compositions"] += 1
        composed = composer.compose_many(outer, argument)
        statistics["batch_calls"] += 1
        values = [[[R.ctx([]), R.ctx([])] for _ in range(L)] for _ in range(L)]
        for (i, j, part), value in zip(positions, composed):
            values[i][j][part] = value
        matrix_out = [[tuple(value) for value in row] for row in values]
        return matrix_out, composed[-1] if include_jump else argument

    try:
        if d > 1:
            x = R.ctx([0, 1]) % R.g
            argument = x.mul_mod(x, R.g)
            statistics["argument_squares"] += 1
            for index, digit in enumerate(bin(d)[3:]):
                shifted, next_argument = compose(product, argument, True, index == 0)
                product = multiply(shifted, product)
                argument = next_argument
                statistics["doublings"] += 1
                statistics["frobenius_matrix_products"] += 1
                if digit == "1":
                    shifted, _ = compose(A, argument, False, False)
                    product = multiply(shifted, product)
                    argument = argument.mul_mod(argument, R.g)
                    statistics["argument_squares"] += 1
                    statistics["appends"] += 1
                    statistics["frobenius_matrix_products"] += 1
            if argument != x:
                raise ArithmeticError("paid Frobenius orbit does not close at d")
    finally:
        if composer is not None:
            composer.close()
            statistics["composer_closed"] = True
            statistics["native"] = dict(composer.statistics_snapshot(),
                                         created=True, closed=True, skip_reason=None)
        else:
            statistics["native"] = dict(
                created=False, closed=True, calls=0,
                skip_reason="degree_one" if d == 1 else "first_inflate_only",
            )

    # The Frobenius cycle is complete. Keep its whole matrix whenever a true
    # N/d power remains; only that power's final multiplication may contract.
    power_base = product
    for digit in bin(e)[3:]:
        product = multiply(product, product)
        statistics["subfield_matrix_products"] += 1
        if digit == "1":
            product = multiply(power_base, product)
            statistics["subfield_matrix_products"] += 1
    if products_done != expected_products:
        raise ArithmeticError("binary product schedule ended early")
    if final_trace is None:
        if product is None:
            raise ArithmeticError("missing both final matrix and complete trace")
        final_trace = R.zero
        for j in range(L):
            final_trace = R.add(final_trace, product[j][j])
    statistics["packed"] = multiplier.statistics_snapshot()
    statistics["parent_xi_constructed"] = parent.xi is not None
    # Releasing these query-owned native and input references is part of the
    # enclosing binary stage. Debug deliberately retains only specified data.
    composer = multiplier = power_base = A = None
    return final_trace, product, statistics, final_pair


def query(acquired, method, N, L, f, a, mode="base", keep_debug=False):
    """Read one actual-N parameter, with no candidate/reference global patching."""
    started = perf_counter()
    _mode(mode)
    if type(keep_debug) is not bool:
        raise ValueError("keep_debug must be a boolean")
    if method not in ("ah", "ordinary") or not isinstance(acquired, dict):
        raise ValueError("invalid acquisition or method")
    if type(N) is not int or N < 1 or type(L) is not int or L < 2:
        raise ValueError("expected N>=1 and L>=2")
    if type(f) is not int or f <= 0 or f.bit_length() != N + 1:
        raise ValueError("expected an explicit monic binary degree-N modulus")
    if type(a) is not int or a <= 0 or a >= 1 << N:
        raise ValueError("expected a nonzero degree-N parameter")
    if acquired.get("N") != N or acquired.get("L") != L:
        raise ValueError("acquisition metadata must use the actual N and L")
    P = acquired.get("P", acquired.get("output_P"))
    if type(P) is not int or P < 1:
        raise ValueError("acquisition must specify its actual output precision")
    Q = acquired.get("Q") if method == "ah" else P
    if type(Q) is not int or Q != (P + L - 1 if method == "ah" else P):
        raise ValueError("AH row division requires the complete registered guard")
    fa, orbit = common.parameter_polynomial(a, f)
    presentation = ForwardPresentation if mode in ("lazy_tau", "combined") else common.Presentation
    parent = presentation(fa, Q)
    prepared = perf_counter()
    B = common.specialize(acquired, method, parent, L, P)
    specialized = perf_counter()
    trace, product, statistics, pair = ordered_trace(
        parent, B, acquired["gamma"], N, L, P, mode, keep_debug)
    value = recover_trace(parent.ring(P), trace, L)
    row = dict(implementation=IMPLEMENTATION, mode=mode, a=a, value=value,
               degree=len(orbit), fa=fa, P=P, Q=Q, operations=statistics)
    if keep_debug:
        row["_debug"] = dict(
            parent=parent, B=B, trace=trace, final_matrix=product,
            terminal_U=pair[0] if pair is not None else None,
            terminal_V=pair[1] if pair is not None else None,
        )
    parent = B = trace = product = pair = orbit = None
    ended = perf_counter()
    row["stages"] = dict(parameter=prepared-started,
                         specialization=specialized-prepared,
                         binary_product=ended-specialized, total=ended-started)
    return row
