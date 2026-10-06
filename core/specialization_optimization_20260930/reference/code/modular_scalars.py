"""Finite-modulus scalar acquisition for the binary arithmetic construction.

Public API: scalar_data(P) -> (gamma, z, phi, b), with the same lengths and
residues as reference/code/vendor/check_arithmetic.py: len(phi) == P and
len(b) == 2*P.  Every returned entry is a Python int in [0, 2**P).

This implements Appendix A, "Obtaining the scalars", with P_work=5*P+16.
All nonunit divisions explicitly consume precision; no Fraction, unrestricted
integer powering, optional native dependency, or cross-call cache is used.
The Newton lift uses the odd derivative of G(z)=sum_r
2**(2**r-r-1)*z**(2**r-1), so a correct m-bit root lifts to min(2*m,P_work)
bits.  SCALAR_OPTIMIZATION.md gives the equivalence and precision argument.
"""


def _hensel_value_derivative(z: int, bits: int) -> tuple[int, int]:
    """Evaluate the finite, exact reduction of G and G' modulo 2**bits."""
    mask = (1 << bits) - 1
    value = (1 + z) & mask
    derivative = 1
    # power=z**(2**(r-1)-1), slope=its derivative, before the next step.
    power = z & mask
    slope = 1
    degree, r = 4, 2
    while degree - r - 1 < bits:
        square = (power * power) & mask
        slope = (square + ((2 * power * slope) & mask) * z) & mask
        power = (square * z) & mask
        valuation = degree - r - 1
        value = (value + (power << valuation)) & mask
        derivative = (derivative + (slope << valuation)) & mask
        degree <<= 1
        r += 1
    return value, derivative


def _hensel_root(work_bits: int) -> int:
    """The unique odd root of G modulo 2**work_bits; no persistent cache."""
    if work_bits < 1:
        raise ValueError("working precision must be positive")
    z, bits = 1, 1
    while bits < work_bits:
        target = min(2 * bits, work_bits)
        modulus = 1 << target
        value, derivative = _hensel_value_derivative(z, target)
        # The derivative is odd because every term after z has even coefficient.
        z = (z - value * pow(derivative, -1, modulus)) & (modulus - 1)
        bits = target
    return z


def _artin_hasse_coefficients(p: int, work_bits: int) -> tuple[list[int], int]:
    """Return b_0,...,b_(2P-1) modulo 2**P and final guaranteed precision.

    After step k, the newest coefficient is known modulo
    2**(work_bits-v2(k!)); earlier coefficients have at least those digits.
    """
    coefficients = [1]
    bits = work_bits
    for k in range(1, 2 * p):
        total = 0
        offset = 1
        while offset <= k:
            total += coefficients[k - offset]
            offset <<= 1
        total &= (1 << bits) - 1
        loss = (k & -k).bit_length() - 1
        divisor = 1 << loss
        if total & (divisor - 1):
            raise ArithmeticError("Artin--Hasse numerator is not divisible by its 2-part")
        bits -= loss
        if bits < p:
            raise ArithmeticError("insufficient guard precision in Artin--Hasse recurrence")
        modulus = 1 << bits
        odd_part = k >> loss
        divided = total >> loss
        if odd_part != 1:
            divided *= pow(odd_part, -1, modulus)
        coefficients.append(divided & (modulus - 1))
    mask = (1 << p) - 1
    return [coefficient & mask for coefficient in coefficients], bits


def _potential_coefficients(p: int, z: int) -> list[int]:
    """Actual phi(y)=y+sum a_ell*y**(2**ell), to degree <P modulo 2**P."""
    phi = [0] * p
    if p > 1:
        phi[1] = 1
    if p <= 2:
        return phi
    modulus = 1 << p
    mask = modulus - 1
    # All requested P+D are <2P.  Keep the common numerator at 2P bits.
    numerator_mask = (1 << (2 * p)) - 1
    gamma_power = (2 * z) & numerator_mask
    numerator = gamma_power  # n_0=gamma; n_ell=2*n_(ell-1)+gamma**(2**ell).
    inverse_z_power = pow(z & mask, -1, modulus)
    degree = 2
    while degree < p:
        gamma_power = (gamma_power * gamma_power) & numerator_mask
        numerator = ((numerator << 1) + gamma_power) & numerator_mask
        inverse_z_power = (inverse_z_power * inverse_z_power) & mask
        residue = numerator & ((1 << (p + degree)) - 1)
        if residue & ((1 << degree) - 1):
            raise ArithmeticError("correction numerator is not divisible by 2**D")
        phi[degree] = ((residue >> degree) * inverse_z_power) & mask
        degree <<= 1
    return phi


def scalar_data(p: int) -> tuple[int, int, list[int], list[int]]:
    """Acquire every scalar from scratch at target precision P=p>=1.

    There is deliberately no lru_cache: cold-start callers pay the Hensel lift,
    actual correction acquisition, and the entire coefficient recurrence.
    """
    if not isinstance(p, int) or isinstance(p, bool):
        raise TypeError("P must be an integer")
    if p < 1:
        raise ValueError("P must be positive")
    work_bits = 5 * p + 16
    z = _hensel_root(work_bits)
    phi = _potential_coefficients(p, z)
    b, _ = _artin_hasse_coefficients(p, work_bits)
    mask = (1 << p) - 1
    return (2 * z) & mask, z & mask, phi, b
