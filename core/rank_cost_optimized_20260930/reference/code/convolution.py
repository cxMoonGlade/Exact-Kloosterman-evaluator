"""Paid, exact finite-field convolution baselines for the rank-cost experiment.

Frozen representation: integers encode polynomial-basis GF(2^N) elements, with
f including its monic x^N bit.  Default polynomials are 0x13, 0x43, 0x11b,
0x1053, and 0x1100b for N=4,6,8,12,16 respectively.  Each run checks
irreducibility, finds a primitive element, and builds its exponent/logarithm
index.  No pre-existing logarithm or multiplication table is used.

Derivation (finite sums, independent of the manuscript's p-adic reduction):
write a=g^k.  Multiplicative convolution becomes cyclic convolution
    (u*v)[k] = sum_{i=0}^{q-2} u[i] v[(k-i) mod (q-1)].
An L-fold convolution of psi therefore sums psi(x_1)...psi(x_L)
over x_1...x_L=a.  Additivity of the absolute trace makes this the defining
hyper-Kloosterman sum.  The convolution identity is delta_1, i.e. index 0.
Associativity permits binary powering.  The sequential implementation uses
L-1 convolutions with the character (whose entries are +/-1); binary powering
uses general integer multiplication, including delta_1's first product.

All integer arithmetic is Python arbitrary precision.  Dense convolution
uses m^2 terms per call for m=q-1, with O(m) arrays.  The binary version has
popcount(L)+floor(log2(L)) such calls, and remains exponential in N.  These
operation counts do not predict speed: factor bit lengths and +/-1 additions
matter.  The implementation neither imports numpy nor uses FFT rounding.

run(..., instrument=False) times only required construction and requested
lookups.  instrument=True adds separate memory/array observations and its
times MUST NOT be used as normal performance timings.  Import/process time
is measured by the external experiment supervisor.  Algorithm time ends
before diagnostic scans and JSON serialization.  Outputs are never cached
across runs.  A single cold query pays for the entire table; later queries
reuse the logarithm index and the full table.
"""

from __future__ import annotations

from time import perf_counter
import sys
import tracemalloc


FIELD_POLYNOMIALS = {4: 0x13, 6: 0x43, 8: 0x11B, 12: 0x1053, 16: 0x1100B}


def field_mul(a: int, b: int, N: int, f: int) -> int:
    """Polynomial-basis multiplication; same encoding as check_arithmetic.py."""
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & (1 << N):
            a ^= f
    return out


def field_pow(a: int, exponent: int, N: int, f: int) -> int:
    out = 1
    while exponent:
        if exponent & 1:
            out = field_mul(out, a, N, f)
        exponent >>= 1
        if exponent:
            a = field_mul(a, a, N, f)
    return out


def _poly_remainder(a: int, b: int) -> int:
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def _poly_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, _poly_remainder(a, b)
    return a


def irreducible(N: int, f: int) -> bool:
    """Finite-field criterion, paid at each cold start (no polynomial table oracle)."""
    if N < 1 or f.bit_length() != N + 1 or not (f & 1):
        return False
    x = _poly_remainder(2, f)
    h = x
    for degree in range(1, N + 1):
        h = field_mul(h, h, N, f)
        if degree <= N // 2 and _poly_gcd(h ^ x, f) != 1:
            return False
    return h == x


def _prime_factors(n: int) -> list[int]:
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


def _cyclic_index(N: int, f: int) -> tuple[int, list[int], list[int]]:
    q = 1 << N
    m = q - 1
    factors = _prime_factors(m)
    primitive = next(
        a for a in range(1, q)
        if all(field_pow(a, m // r, N, f) != 1 for r in factors)
    )
    exponents = [0] * m
    logarithms = [-1] * q
    a = 1
    for i in range(m):
        exponents[i] = a
        logarithms[a] = i
        a = field_mul(a, primitive, N, f)
    if a != 1:
        raise ArithmeticError("Primitive-element cycle failed to close")
    return primitive, exponents, logarithms


def absolute_trace(a: int, N: int, f: int) -> int:
    trace = 0
    for _ in range(N):
        trace ^= a
        a = field_mul(a, a, N, f)
    if trace not in (0, 1):
        raise ArithmeticError("Absolute trace is outside GF(2)")
    return trace


def cyclic_general(a: list[int], b: list[int]) -> list[int]:
    """One general dense cyclic convolution; m^2 arbitrary-int multiply/adds."""
    m = len(a)
    if len(b) != m:
        raise ValueError("Convolution lengths differ")
    out = [0] * m
    for i, x in enumerate(a):
        k = i
        for y in b:
            out[k] += x * y
            k += 1
            if k == m:
                k = 0
    return out


def cyclic_sign(a: list[int], signs: list[int]) -> list[int]:
    """Dense cyclic convolution with a fixed +/-1 second factor."""
    m = len(a)
    if len(signs) != m:
        raise ValueError("Convolution lengths differ")
    out = [0] * m
    for i, x in enumerate(a):
        k = i
        for sign in signs:
            if sign == 1:
                out[k] += x
            else:
                out[k] -= x
            k += 1
            if k == m:
                k = 0
    return out


def _array_summary(values: list[int]) -> dict:
    return {
        "length": len(values),
        "nonzero_elements": sum(x != 0 for x in values),
        "max_integer_bits": max((abs(x).bit_length() for x in values), default=0),
        "container_bytes": sys.getsizeof(values),
    }


def _retained_python_bytes(arrays: tuple[list[int], ...]) -> int:
    """Actual retained list/int objects, counting aliased Python objects once."""
    seen = set()
    total = 0
    for values in arrays:
        if id(values) not in seen:
            seen.add(id(values))
            total += sys.getsizeof(values)
        for value in values:
            if id(value) not in seen:
                seen.add(id(value))
                total += sys.getsizeof(value)
    return total


def run(
    N: int,
    L: int,
    f: int,
    params: list[int],
    instrument: bool = False,
    method: str = "binary",
) -> dict:
    """Return exact requested integers and disjoint construction/query timings.

    preprocess_seconds includes field validation, index construction, character,
    and the full convolution table.  cold_seconds is preprocessing plus the
    first query; compute_seconds is preprocessing plus ALL requested queries
    and small stage-boundary overhead.  Each requested parameter is a nonzero
    polynomial-basis element.  Metadata/memory observations occur only when
    instrument=True; those runs are explicitly labelled timing_eligible=False.
    """
    if method not in ("sequential", "binary"):
        raise ValueError("method must be 'sequential' or 'binary'")
    if L < 1 or N < 1:
        raise ValueError("N and L must be positive")
    if not params or any(not 0 < a < 1 << N for a in params):
        raise ValueError("At least one nonzero field parameter is required")
    own_trace = instrument and not tracemalloc.is_tracing()
    if own_trace:
        tracemalloc.start()
    snapshots = []
    start = perf_counter()
    if not irreducible(N, f):
        raise ValueError("The supplied modulus is not irreducible of degree N")
    primitive, exponents, logarithms = _cyclic_index(N, f)
    after_index = perf_counter()
    psi = [1 - 2 * absolute_trace(a, N, f) for a in exponents]
    if instrument:
        snapshots.append({
            "phase": "index_and_character",
            "arrays": {name: _array_summary(a) for name, a in
                       (("exponents", exponents), ("logarithms", logarithms), ("character", psi))},
        })
    del exponents
    after_character = perf_counter()
    m = (1 << N) - 1
    convolution_calls = 0
    maximum_observed_value_bits = 1
    if method == "sequential":
        table = psi
        for _ in range(1, L):
            out = cyclic_sign(table, psi)
            if instrument:
                convolution_calls += 1
                snapshots.append({
                    "phase": f"convolution_{convolution_calls}",
                    "arrays": {name: _array_summary(a) for name, a in
                               (("left", table), ("right_character", psi), ("output", out))},
                })
                maximum_observed_value_bits = max(
                    maximum_observed_value_bits, snapshots[-1]["arrays"]["output"]["max_integer_bits"])
            table = out
        del psi
    else:
        base = psi
        del psi
        table = [1] + [0] * (m - 1)
        exponent = L
        while exponent:
            if exponent & 1:
                out = cyclic_general(table, base)
                if instrument:
                    convolution_calls += 1
                    snapshots.append({
                        "phase": f"multiply_{convolution_calls}",
                        "arrays": {name: _array_summary(a) for name, a in
                                   (("accumulator", table), ("base", base), ("output", out))},
                    })
                    maximum_observed_value_bits = max(
                        maximum_observed_value_bits, snapshots[-1]["arrays"]["output"]["max_integer_bits"])
                table = out
            exponent >>= 1
            if exponent:
                out = cyclic_general(base, base)
                if instrument:
                    convolution_calls += 1
                    snapshots.append({
                        "phase": f"square_{convolution_calls}",
                        "arrays": {name: _array_summary(a) for name, a in
                                   (("accumulator", table), ("base", base), ("output", out))},
                    })
                    maximum_observed_value_bits = max(
                        maximum_observed_value_bits, snapshots[-1]["arrays"]["output"]["max_integer_bits"])
                base = out
        del base
    # Avoid keeping an otherwise unused alias to the final convolution output.
    if method == "binary" or L > 1:
        del out
    after_table = perf_counter()
    outputs = []
    query_seconds = []
    for a in params:
        query_start = perf_counter()
        value = table[logarithms[a]]
        query_end = perf_counter()
        outputs.append(value)
        query_seconds.append(query_end - query_start)
    end = perf_counter()
    stages = {
        "field_and_index": after_index - start,
        "character": after_character - after_index,
        "convolution": after_table - after_character,
        "lookup": sum(query_seconds),
    }
    preprocessing = after_table - start
    result = {
        "method": f"convolution_{method}",
        "N": N, "L": L, "f": f, "f_hex": hex(f), "params": params,
        "outputs": outputs, "primitive_element": primitive,
        "query_count": len(params), "full_table_entries": m,
        "preprocess_seconds": preprocessing,
        "query_seconds": query_seconds,
        "cold_seconds": preprocessing + query_seconds[0],
        "compute_seconds": end - start,
        "stages_seconds": stages,
        "stage_gap_seconds": end - start - sum(stages.values()),
        "timing_eligible": not instrument,
        "instrumented": instrument,
        "table_is_all_parameters": True,
        "reusable_objects": ["logarithm_index", "full_convolution_table"],
    }
    if instrument:
        sign_calls = convolution_calls if method == "sequential" else 0
        current_bytes, peak_bytes = tracemalloc.get_traced_memory()
        # The allocation peak is sampled BEFORE potentially large identity sets
        # used by the final exact retained-object census.
        if own_trace:
            tracemalloc.stop()
        result["instrumentation"] = {
            "tracemalloc_current_before_census_bytes": current_bytes,
            "tracemalloc_peak_before_census_bytes": peak_bytes,
            "tracemalloc_scope": "this_run" if own_trace else "external_tracer",
            "retained_array_slots": len(logarithms) + len(table),
            "retained_arrays_python_bytes_deduplicated": _retained_python_bytes((logarithms, table)),
            "retained_arrays": {
                "logarithm_index": _array_summary(logarithms),
                "full_convolution_table": _array_summary(table),
            },
            "maximum_observed_value_bits": maximum_observed_value_bits,
            "convolution_calls": convolution_calls,
            "sign_convolution_calls": sign_calls,
            "general_convolution_calls": convolution_calls - sign_calls,
            "integer_additions_derived_from_executed_dense_loops": convolution_calls * m * m,
            "integer_multiplications_derived_from_executed_dense_loops": (convolution_calls - sign_calls) * m * m,
            "array_snapshots": snapshots,
            "note": "Memory includes instrumented algorithm and small observation metadata; final deduplicating census excluded from allocation peak. Snapshot containers are measured; no theoretical index count is presented as RAM.",
        }
    return result


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=4)
    parser.add_argument("--L", type=int, default=4)
    parser.add_argument("--f", type=lambda s: int(s, 0))
    parser.add_argument("--params", type=int, nargs="+", default=[1, 2, 3, 7])
    parser.add_argument("--method", choices=["sequential", "binary"], default="binary")
    parser.add_argument("--instrument", action="store_true")
    args = parser.parse_args()
    modulus = args.f if args.f is not None else FIELD_POLYNOMIALS[args.N]
    print(json.dumps(run(args.N, args.L, modulus, args.params, args.instrument, args.method), indent=2))
