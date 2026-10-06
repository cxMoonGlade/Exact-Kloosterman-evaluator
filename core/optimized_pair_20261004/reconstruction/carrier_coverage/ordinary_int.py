"""D10 Python-integer ordinary Bessel acquisition.

Explicit source copy of frozen code/ordinary.py (SHA f84437fb...da880).
Only the MPZ scalar constructors and implementation identity change. The
FLINT initial jet, precision schedule, recurrence and checks remain intact.

The initial jet uses an exact integer polynomial product tree.  Subsequent
coefficients use the triangular differential recurrence, with the integral
sector cutoff from ORDINARY_SECTOR_TRANSFER and the odd-rank Taylor bridge.
No Artin--Hasse acquisition or old ordinary evaluator is imported.

``acquire`` returns B = A diag(2**(-j)), not the row-normalized old matrix.
The common flag code applies diag(2**j) on the left implicitly.  Odd ranks use
two integer components for Z_2[i]; the common coefficient Frobenius fixes i.
There is no scientific computation on import.
"""

from time import perf_counter

from flint import fmpz, fmpz_mod_poly_ctx, fmpz_poly


IMPLEMENTATION = "ordinary_sector_binary_split_python_int_d10_v1"


def _positive_int(value, name):
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise ValueError(f"{name} must be a positive integer")


def precision_schedule(N, L, P=None):
    """The certified sector cutoff, including every column-division guard."""
    _positive_int(N, "N")
    _positive_int(L, "L")
    if L < 2:
        raise ValueError("Ordinary matrix acquisition requires L >= 2")
    required = N * (L - 1) + 2
    P = required if P is None else P
    _positive_int(P, "P")
    if P < required:
        raise ValueError("P is below the registered exact-integer precision")
    e = P + L
    R = 2 * ((e + L - 1) // L)
    M = e + (2 * L - 1) * R + 2
    loss = (2 * L - 1) * ((R - 1) - (R - 1).bit_count())
    return {
        "N": N, "L": L, "P": P, "minimum_output_P": required,
        "e": e, "R": R, "M": M,
        "initializer_scaled_bits": M + 2 * (L - 1),
        "initializer_terms": 4 * (M + 2 * (L - 1)) + 1,
        "recurrence_loss": loss, "recurrence_final_bits": M - loss,
        "maximum_column_division": L - 1,
        "gaussian_components": 2 if L & 1 else 1,
        "old_matrix_scale_B": 0, "old_cycle_scale_NB": 0,
    }


def _integer_sum_jet(length, bits):
    """Return G_T(4v) mod (v**length, 2**bits) and product-tree counters.

Each tree block stores the six nontrivial polynomial entries and scalar d
of [[a,b,0],[c,d,0],[s,t,denominator]].  The lower-right entry is a scalar;
the two identically zero entries are never materialized.  Products are exact
over Z[v]/v**length, and the full rational denominator is removed only once.
"""
    _positive_int(length, "length")
    _positive_int(bits, "bits")
    K = 4 * bits + 1
    counters = {"integer_polynomial_products": 0,
                "product_tree_nodes": 0, "product_tree_depth": 0}

    def product(a, b):
        if not a or not b:
            return fmpz_poly([])
        counters["integer_polynomial_products"] += 1
        return a.mul_low(b, length)

    def tree(lo, hi, depth):
        counters["product_tree_nodes"] += 1
        counters["product_tree_depth"] = max(counters["product_tree_depth"], depth)
        if hi - lo == 1:
            j = lo
            denominator = fmpz((2 * j + 2) * (2 * j + 1))
            x = fmpz_poly([-j, 4])
            previous_x = fmpz_poly([1 - j, 4])
            a = (-4 * (2 * j + 1) * x).truncate(length)
            b = (-4 * product(x, previous_x)).truncate(length)
            d = fmpz_poly([denominator])
            zero = fmpz_poly([])
            return a, b, d, zero, d, zero, denominator
        middle = (lo + hi) // 2
        left = tree(lo, middle, depth + 1)
        right = tree(middle, hi, depth + 1)
        a, b, c, d, s, t, denominator_l = left
        A, B, C, D, S, T, denominator_r = right
        # Column states evolve by C_(hi-1) ... C_lo: right times left.
        return (
            product(A, a) + product(B, c),
            product(A, b) + product(B, d),
            product(C, a) + product(D, c),
            product(C, b) + product(D, d),
            product(S, a) + product(T, c) + denominator_r * s,
            product(S, b) + product(T, d) + denominator_r * t,
            denominator_r * denominator_l,
        )

    block = tree(0, K, 0)
    numerator, denominator = block[4], int(block[6])
    # Product_j (2j+2)(2j+1) = (2K)!; Legendre counts its exact 2-part.
    twos = 2 * K - (2 * K).bit_count()
    if denominator & ((1 << twos) - 1) or not (denominator >> twos) & 1:
        raise ArithmeticError("Unexpected binary-splitting denominator valuation")
    modulus = 1 << bits
    odd_inverse = pow(denominator >> twos, -1, modulus)
    divisor_mask = (1 << twos) - 1
    result = []
    for k in range(length):
        value = int(numerator[k])
        if value & divisor_mask:
            raise ArithmeticError("Nonintegral scaled ordinary splitting sum")
        result.append(((value >> twos) * odd_inverse) % modulus)
    counters.update(terms=K, jet_length=length, scaled_bits=bits,
                    denominator_twos=twos,
                    denominator_bits=denominator.bit_length())
    return result, counters


def initial_jet(L, bits, *, return_metadata=False):
    """Return [u**k] G_T(u)**L, k < L, as canonical bits-bit residues.

This function constructs the full off-diagonal initial data.  It exposes a
small independent testing surface without constructing a field or a query.
"""
    _positive_int(L, "L")
    _positive_int(bits, "bits")
    scaled_bits = bits + 2 * (L - 1)
    scaled, metadata = _integer_sum_jet(L, scaled_bits)
    context = fmpz_mod_poly_ctx(1 << scaled_bits)
    powered = context(scaled).pow_trunc(L, L)
    mask = (1 << bits) - 1
    result = []
    for k in range(L):
        value = int(powered[k])
        shift = 2 * k
        if value & ((1 << shift) - 1):
            raise ArithmeticError("Ordinary initial jet is not integrally rescalable")
        result.append((value >> shift) & mask)
    metadata.update(output_bits=bits, power=L)
    return (result, metadata) if return_metadata else result


def _divide_integral(num, denominator, input_bits):
    """Exact modular division of the complete numerator; lose only v2(d)."""
    even = (denominator & -denominator).bit_length() - 1
    bits = input_bits - even
    if bits < 1:
        raise ArithmeticError("Ordinary recurrence exhausted its precision")
    if num & ((1 << even) - 1):
        raise ArithmeticError("Nonintegral ordinary recurrence numerator")
    num >>= even
    odd = denominator >> even
    if odd == 1:
        return num
    # A small odd division avoids multiplying two full-precision residues.
    correction = (-int(num % odd) * pow(pow(2, bits, odd), -1, odd)) % odd
    return (num + (int(correction) << bits)) // odd


def _coefficient_table(jet, precision):
    """Phi_m in signed reversal coordinates, with full certified e bits."""
    L, R, bits, e = (precision[key] for key in ("L", "R", "M", "e"))
    gaussian = bool(L & 1)
    keep_mask = (1 << e) - 1
    initial_mask = (1 << bits) - 1
    previous_real = [
        (int(jet[j - i]) * (-1 if (i + j) & 1 else 1)
         * (1 << (L - 1 - j))) & initial_mask
        if j >= i else 0 for i in range(L) for j in range(L)
    ]
    previous_imag = [0 for _ in range(L * L)]
    older_real = [0 for _ in range(L * L)]
    older_imag = [0 for _ in range(L * L)]
    real = [[value & keep_mask] for value in previous_real]
    imag = [[0] for _ in range(L * L)]
    c = 1 << L
    divisions = 0
    for m in range(1, R):
        rhs_real = [0 for _ in range(L * L)]
        rhs_imag = [0 for _ in range(L * L)]
        for j in range(L):
            target = (L - 1) * L + j
            if gaussian:
                rhs_real[target] -= 2 * c * older_imag[j]
                rhs_imag[target] += 2 * c * older_real[j]
            else:
                rhs_real[target] += 2 * c * older_real[j]
        for i in range(L):
            source, target = i * L + L - 1, i * L
            if gaussian:
                rhs_real[target] += c * previous_imag[source]
                rhs_imag[target] -= c * previous_real[source]
            else:
                rhs_real[target] -= c * previous_real[source]
        new_real = [0 for _ in range(L * L)]
        new_imag = [0 for _ in range(L * L)]
        even = (m & -m).bit_length() - 1
        for i in range(L - 1, -1, -1):
            for j in range(L):
                index = i * L + j
                input_bits = bits - (L - i + j - 1) * even
                mask = (1 << input_bits) - 1
                value = rhs_real[index]
                if j:
                    value -= new_real[index - 1]
                if i + 1 < L:
                    value += 2 * new_real[index + L]
                new_real[index] = _divide_integral(value & mask, m, input_bits)
                divisions += 1
                if gaussian:
                    value = rhs_imag[index]
                    if j:
                        value -= new_imag[index - 1]
                    if i + 1 < L:
                        value += 2 * new_imag[index + L]
                    new_imag[index] = _divide_integral(value & mask, m, input_bits)
                    divisions += 1
        bits -= (2 * L - 1) * even
        if bits < e:
            raise ArithmeticError("Ordinary recurrence fell below acquisition precision")
        common_mask = (1 << bits) - 1
        new_real = [value & common_mask for value in new_real]
        new_imag = [value & common_mask for value in new_imag]
        for index in range(L * L):
            real[index].append(new_real[index] & keep_mask)
            imag[index].append(new_imag[index] & keep_mask)
        older_real, previous_real = previous_real, new_real
        older_imag, previous_imag = previous_imag, new_imag
    if bits != precision["recurrence_final_bits"]:
        raise ArithmeticError("Ordinary recurrence precision accounting mismatch")
    return real, imag, {"recurrence_divisions": divisions,
                        "retained_coefficient_matrices": R,
                        "final_working_bits": bits}


def acquire(N, L, P=None):
    """Acquire nested B_coefficients[i][j] = (real_list, imag_list).

Every list has length R (degrees 0 through R-1), modulo 2**output_P.
The returned gamma is 2 and B = A diag(2**(-j)); B mod 2 is I.
``P`` may increase the standard output precision, e.g. for P versus P+8
checks.  Acquisition is independent of the field polynomial and parameter.
"""
    start = perf_counter()
    precision = precision_schedule(N, L, P)
    stage = perf_counter()
    jet, initializer_counters = initial_jet(L, precision["M"], return_metadata=True)
    stages = {"ordinary_initializer": perf_counter() - stage}
    stage = perf_counter()
    phi_real, phi_imag, recurrence_counters = _coefficient_table(jet, precision)
    stages["ordinary_coefficient_recurrence"] = perf_counter() - stage
    stage = perf_counter()
    e_mask = (1 << precision["e"]) - 1
    output_mask = (1 << precision["P"]) - 1
    result = []
    checked = 0
    for i in range(L):
        row = []
        for j in range(L):
            index = (L - 1 - i) * L + L - 1 - j
            sign = -1 if (i + j) & 1 else 1
            components = []
            for series in (phi_real[index], phi_imag[index]):
                output = []
                for m, value in enumerate(series):
                    value = (sign * value) & e_mask
                    # Certified sector bound before dividing column j.
                    lower = L * (m // 2) + j if m % 2 == 0 else L * (m // 2 + 1)
                    visible_lower = min(lower, precision["e"])
                    if value & ((1 << visible_lower) - 1):
                        raise ArithmeticError("Ordinary sector coefficient filter failed")
                    normalized = (value >> j) & output_mask
                    output.append(int(normalized))
                    checked += 1
                components.append(output)
            real, imag = components
            if (real[0] & 1) != int(i == j) or any(value & 1 for value in real[1:]):
                raise ArithmeticError("Ordinary B is not the constant identity modulo 2")
            if any(value & 1 for value in imag):
                raise ArithmeticError("Ordinary Gaussian B has a nonzero imaginary mod-2 part")
            row.append((real, imag))
            phi_real[index] = phi_imag[index] = None
        result.append(row)
    stages["ordinary_column_normalization"] = perf_counter() - stage
    counters = dict(initializer_counters)
    counters.update(recurrence_counters, checked_sector_coefficients=checked,
                    output_coefficient_slots=2 * L * L * precision["R"])
    return {
        "implementation": IMPLEMENTATION, "N": N, "L": L,
        "gamma": 2, "output_P": precision["P"],
        "B_coefficients": result, "precision": precision,
        "stages": stages, "counters": counters,
        "seconds": perf_counter() - start,
        "coefficient_layout": "nested row-major L x L; each entry (real, imag), ascending degrees",
        "matrix_convention": "B=A*diag(2**(-j)); flag matrix diag(2**j)*B",
        "gaussian": bool(L & 1),
    }


def validate_components():
    """Bounded exact checks, called explicitly by the registered validator.

The jet oracle expands exp(-2y) exp(2y**2) as rational finite sums and
multiplies falling factorials literally.  The ODE oracle uses complete
Gaussian matrix products, independently of the triangular solver order.
"""
    from fractions import Fraction
    from math import factorial

    def polynomial_product(left, right, length):
        out = [Fraction(0) for _ in range(length)]
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                if i + j < length:
                    out[i + j] += a * b
        return out

    jet_cases = [(2, 3), (3, 3), (4, 2)]
    jet_coefficients = 0
    for L, bits in jet_cases:
        T, K = bits + 2 * (L - 1), 4 * (bits + 2 * (L - 1)) + 1
        falling = [Fraction(1)] + [Fraction(0) for _ in range(L - 1)]
        total = [Fraction(0) for _ in range(L)]
        for j in range(K):
            if j:
                falling = polynomial_product(falling, [-(j - 1), 4], L)
            degree = 2 * j
            beta = sum((Fraction((-2) ** (degree - 2 * k) * 2 ** k,
                                 factorial(degree - 2 * k) * factorial(k))
                        for k in range(j + 1)), Fraction(0))
            factor = beta / ((-2) ** j)
            total = [a + factor * b for a, b in zip(total, falling)]
        powered = [Fraction(1)] + [Fraction(0) for _ in range(L - 1)]
        for _ in range(L):
            powered = polynomial_product(powered, total, L)
        expected = []
        modulus = 1 << bits
        for k, value in enumerate(powered):
            value /= 4 ** k
            if value.denominator % 2 == 0:
                raise AssertionError("Literal ordinary jet oracle is not 2-integral")
            expected.append((value.numerator * pow(value.denominator, -1, modulus)) % modulus)
        actual = initial_jet(L, bits)
        if actual != expected:
            raise AssertionError({"ordinary_initial_jet": (L, bits),
                                  "actual": actual, "expected": expected, "tail_bits": T})
        jet_coefficients += L

    def add(a, b):
        return a[0] + b[0], a[1] + b[1]

    def multiply(a, b):
        return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]

    def matrix_product(A, B):
        n = len(A)
        out = [[(0, 0) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    out[i][j] = add(out[i][j], multiply(A[i][k], B[k][j]))
        return out

    recurrence_cases = [(1, 2), (2, 3), (3, 4), (2, 5)]
    residual_components = 0
    filter_components = 0
    for N, L in recurrence_cases:
        precision = precision_schedule(N, L)
        jet = initial_jet(L, precision["M"])
        real, imag, counters = _coefficient_table(jet, precision)
        if counters["final_working_bits"] < precision["e"]:
            raise AssertionError("Insufficient final ordinary coefficient precision")
        zero = [[(0, 0) for _ in range(L)] for _ in range(L)]
        J = [[(int(j == i + 1), 0) for j in range(L)] for i in range(L)]
        C = [[(0, 0) for _ in range(L)] for _ in range(L)]
        C[L - 1][0] = (0, 1 << L) if L & 1 else (1 << L, 0)
        matrices = [
            [[(int(real[i * L + j][m]), int(imag[i * L + j][m]))
              for j in range(L)] for i in range(L)]
            for m in range(precision["R"])
        ]
        mask = (1 << precision["e"]) - 1
        for m, current in enumerate(matrices):
            previous = matrices[m - 1] if m >= 1 else zero
            older = matrices[m - 2] if m >= 2 else zero
            products = [matrix_product(current, J), matrix_product(J, current),
                        matrix_product(C, older), matrix_product(previous, C)]
            for i in range(L):
                for j in range(L):
                    for component in range(2):
                        residual = (m * current[i][j][component]
                                    + products[0][i][j][component]
                                    - 2 * products[1][i][j][component]
                                    - 2 * products[2][i][j][component]
                                    + products[3][i][j][component])
                        if residual & mask:
                            raise AssertionError({"ordinary_ODE": (N, L, m, i, j, component),
                                                  "residual": residual & mask})
                        residual_components += 1
        acquired = acquire(N, L)
        output_mask = (1 << precision["P"]) - 1
        for i in range(L):
            for j in range(L):
                for component in range(2):
                    for m, value in enumerate(acquired["B_coefficients"][i][j][component]):
                        raw = matrices[m][L - 1 - i][L - 1 - j][component]
                        raw = ((-1) ** (i + j) * raw) & mask
                        if raw & ((1 << j) - 1) or value != ((raw >> j) & output_mask):
                            raise AssertionError("Ordinary signed reversal/column normalization mismatch")
                        filter_components += 1
    return {"status": "PASS", "initial_jet_cases": jet_cases,
            "initial_jet_coefficients": jet_coefficients,
            "recurrence_cases": recurrence_cases,
            "complete_ODE_residual_components": residual_components,
            "independent_column_normalizations": filter_components,
            "oracle": "literal rational exponential/falling sum and full Gaussian matrix products"}
