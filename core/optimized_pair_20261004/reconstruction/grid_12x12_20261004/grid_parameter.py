"""12x12 execution envelope: frozen D14 constant-coordinate algorithm.

These are internal primitives with an already certified irreducible binary
field as a precondition. The registered execution boundary calls
common.validate_input and pays for it. checked_input only checks types,
shapes, supported degree, and capacities; it is not an irreducibility test.
No mathematical operation occurs on import.
"""
from struct import calcsize
import sys
from time import perf_counter

import flint
from flint import fmpz_mod_poly_ctx


IMPLEMENTATION = "constant_coordinate_minpoly_flint_v1"
MAXIMUM_N = 128
_STAGES = ("check", "encode", "orbit", "powers", "recurrence", "output", "cleanup")
_COUNTERS = ("orbit_squares", "orbit_comparisons", "power_products",
             "mul_mod_calls", "sequence_getters", "output_coefficients", "ctx_minpoly_calls")


def checked_input(N, f, a):
    """Reject unsupported dimensions before shifts, arrays, or native factories.

    The caller must independently certify f's irreducibility at the execution
    boundary. No bool flag or purported certification record is accepted here.
    """
    if type(N) is not int:
        raise TypeError("N must be a Python integer, not bool")
    if not 1 <= N <= MAXIMUM_N:
        raise ValueError("Grid supports 1 <= N <= 128; larger degrees require a new review")
    if type(f) is not int or type(a) is not int:
        raise TypeError("f and a must be Python integers, not bool")
    if f <= 0 or f.bit_length() != N + 1:
        raise ValueError("expected a monic degree-N binary polynomial")
    if a <= 0 or a.bit_length() > N:
        raise ValueError("expected a nonzero parameter in the original degree-N basis")
    if tuple(calcsize(fmt) for fmt in ("l", "L", "N", "P", "i")) != (8, 8, 8, 8, 4):
        raise RuntimeError("D14 requires the reviewed LP64 platform with a 32-bit C int")

    # All arithmetic here is Python integer arithmetic, before any C narrowing.
    signed_limit = min(sys.maxsize, (1 << 63) - 1)
    size_limit = (1 << 64) - 1
    int_limit = (1 << 31) - 1
    sequence_bound = 2 * N
    if sequence_bound > int_limit:
        raise OverflowError("sequence length exceeds the native index bound")
    # Input lists, the sequence, the returned polynomial, mulmod's combined
    # product/quotient workspace, and conservative fit_length doubling.
    slots = (N + 1, N, sequence_bound, sequence_bound + 1,
             2 * N - 1, 3 * N - 2, 2 * (sequence_bound + 1))
    for length in slots:
        if length < 1 or length > signed_limit:
            raise OverflowError("native polynomial length exceeds slong/Py_ssize_t")
        if length * 8 > min(signed_limit, size_limit):
            raise OverflowError("coefficient or list storage exceeds the capacity contract")
    for byte_length in ((N + 8) // 8, (N + 7) // 8):
        if byte_length < 1 or byte_length > signed_limit:
            raise OverflowError("input byte buffer exceeds the capacity contract")
    # With modulus 2, FLINT 3.6's BM cutoff is max(200, 530-22*1)=508.
    # N <= 128 gives sequence_bound <= 256, so HGCD remains unreachable.
    return None


def _binary_coefficients(value, length):
    """Little-endian bits without successively shifting a long Python integer."""
    data = value.to_bytes((length + 7) // 8, "little")
    coefficients = []
    for byte in data:
        for bit in range(8):
            if len(coefficients) == length:
                return coefficients
            coefficients.append((byte >> bit) & 1)
    return coefficients


def parameter_polynomial(N, f, a):
    """Return (fa, actual_degree, detail) for an already certified field input.

    fa is ordered from constant coefficient to leading coefficient. The
    method pays for d actual squarings, 2*d-1 power products, and one length
    2*d minimal-recurrence solve, including d=1. It does not return an orbit.
    """
    started = phase_started = perf_counter()
    phase = "check"
    seconds = {name: 0.0 for name in _STAGES}
    seconds["total"] = 0.0
    counts = {name: 0 for name in _COUNTERS}
    R = F = A = current = power = recurrence = None
    f_coefficients = a_coefficients = sequence = None
    fa = detail = None
    d = 0

    def next_phase(name):
        nonlocal phase, phase_started
        now = perf_counter()
        seconds[phase] += now - phase_started
        phase, phase_started = name, now

    try:
        checked_input(N, f, a)
        next_phase("encode")
        if (flint.__version__, flint.__FLINT_VERSION__, flint.ctx.threads) != ("0.9.0", "3.6.0", 1):
            raise RuntimeError("D14 requires the reviewed single-thread python-flint/FLINT versions")
        f_coefficients = _binary_coefficients(f, N + 1)
        a_coefficients = _binary_coefficients(a, N)
        R = fmpz_mod_poly_ctx(2)
        F, A = R(f_coefficients), R(a_coefficients)
        current = A
        power = R([1])

        next_phase("orbit")
        for step in range(1, N + 1):
            current = current.mul_mod(current, F)
            counts["orbit_squares"] += 1
            counts["mul_mod_calls"] += 1
            returned = current == A
            counts["orbit_comparisons"] += 1
            if returned:
                d = step
                break
        if d < 1 or d > N or N % d:
            raise ArithmeticError("squaring did not give a valid first-return period dividing N")

        next_phase("powers")
        sequence = []
        for j in range(2 * d):
            coefficient = int(power[0])
            counts["sequence_getters"] += 1
            if coefficient not in (0, 1):
                raise ArithmeticError("constant-coordinate sequence is not binary")
            sequence.append(coefficient)
            if j + 1 < 2 * d:
                power = power.mul_mod(A, F)
                counts["power_products"] += 1
                counts["mul_mod_calls"] += 1

        next_phase("recurrence")
        recurrence = R.minpoly(sequence)
        counts["ctx_minpoly_calls"] += 1

        next_phase("output")
        if recurrence.degree() != d:
            raise ArithmeticError("minimal recurrence degree differs from the actual orbit degree")
        fa = []
        for j in range(d + 1):
            coefficient = int(recurrence[j])
            counts["output_coefficients"] += 1
            if coefficient not in (0, 1):
                raise ArithmeticError("minimal recurrence coefficient is not binary")
            fa.append(coefficient)
        if fa[-1] != 1:
            raise ArithmeticError("minimal recurrence is not monic")
        detail = dict(method=IMPLEMENTATION, N=N, degree=d, sequence_length=len(sequence),
                      counts=counts, seconds=seconds)
    finally:
        # Closing the active bucket before cleanup also accounts for early
        # rejection/exception time. No native result or context is in detail.
        next_phase("cleanup")
        recurrence = power = current = A = F = R = None
        f_coefficients = a_coefficients = sequence = None
        ended = perf_counter()
        seconds["cleanup"] += ended - phase_started
        seconds["total"] = ended - started
    return fa, d, detail
