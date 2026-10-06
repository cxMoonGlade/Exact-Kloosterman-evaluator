"""Exact Kronecker-packed cyclic convolution, with paid field construction.

This is an accelerated member of the multiplicative-convolution family, not
an independent cohomology/Bessel/deformation algorithm.  No floating-point FFT
or fixed-width integer arithmetic is used.  See ../FAST_CONVOLUTION.md for the
complete finite coefficient/carry bound and task/cost interpretation.

Reference ownership is explicit: ../reference/code/convolution.py supplies
only paid finite-field validation/index/character primitives.  run never uses
its convolution routine, precomputed answer, or cache.  The independent small
literal check is loaded only by --check, never during performance runs.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import struct
import sys
from time import perf_counter
import tracemalloc


REFERENCE_CODE = Path(__file__).resolve().parents[1] / "reference" / "code"
MEMORY_LIMIT_BYTES = 16 * 1024**3


def _load_file(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load frozen reference: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


reference = _load_file("_qcb_fast_convolution_reference", REFERENCE_CODE / "convolution.py")
FIELD_POLYNOMIALS = reference.FIELD_POLYNOMIALS


def memory_precheck(N: int, limit: int = MEMORY_LIMIT_BYTES) -> dict:
    """Strict necessary pointer storage for this implementation, no allocations.

    The retained logarithm index has q slots and answer table q-1 slots.
    Their pointer bytes alone already exceed 16 GiB at N=32 on 64-bit Python,
    and likewise for larger N.  Integer payloads, list headers,
    construction and multiplication temporaries only add to this lower bound.
    This is an implementation feasibility bound, not a universal lower bound.
    """
    if N < 1:
        raise ValueError("N must be positive")
    q = 1 << N
    pointer_size = struct.calcsize("P")
    necessary_bytes = (2 * q - 1) * pointer_size
    return {
        "feasible_by_pointer_lower_bound": necessary_bytes <= limit,
        "required_retained_pointer_bytes_lower_bound": necessary_bytes,
        "pointer_size_bytes": pointer_size,
        "retained_array_slots": 2 * q - 1,
        "memory_limit_bytes": limit,
        "scope": "this full-table Python-list implementation; necessary bound only",
    }


def _pack(values: list[int], offset: int, width: int) -> int:
    """Serialize nonnegative base-2^(8*width) digits without quadratic shifts."""
    raw = bytearray(len(values) * width)
    for i, value in enumerate(values):
        raw[i * width:(i + 1) * width] = (value + offset).to_bytes(width, "little")
    return int.from_bytes(raw, "little")


def cyclic_kronecker(a: list[int], b: list[int], observations: list | None = None) -> list[int]:
    """Exact dense cyclic convolution with one packed arbitrary-int product.

    Let alpha=max(0,-min(a)), beta=max(0,-min(b)), A=a+alpha,
    B=b+beta and C=m*max(A)*max(B).  Every coefficient of the ordinary
    polynomial product A(X)B(X) lies in [0,C].  Choose byte-aligned w with
    2^w>max(C,max(A),max(B)); the latter two bounds also handle zero products.
    Consequently integer multiplication at X=2^w has no carries
    between coefficient blocks.  Fold degree k+m into k, then subtract
    beta*sum(a)+alpha*sum(b)+m*alpha*beta from every cyclic coefficient.
    """
    m = len(a)
    if not m or len(b) != m:
        raise ValueError("Nonempty convolution operands must have equal lengths")
    alpha = max(0, -min(a))
    maximum_a = max(a) + alpha
    sum_a = sum(a)
    if a is b:
        beta, maximum_b, sum_b = alpha, maximum_a, sum_a
    else:
        beta = max(0, -min(b))
        maximum_b = max(b) + beta
        sum_b = sum(b)
    coefficient_bound = m * maximum_a * maximum_b
    digit_bound = max(coefficient_bound, maximum_a, maximum_b)
    width = max(1, (digit_bound.bit_length() + 7) // 8)
    raw_product_bytes = (2 * m - 1) * width
    if raw_product_bytes > MEMORY_LIMIT_BYTES:
        raise MemoryError("Packed product byte buffer alone exceeds configured memory cap")
    correction = beta * sum_a + alpha * sum_b + m * alpha * beta
    packed_a = _pack(a, alpha, width)
    packed_b = packed_a if a is b else _pack(b, beta, width)
    packed_product = packed_a * packed_b
    if observations is not None:
        observation = {
            "length": m,
            "offset_a": alpha,
            "offset_b": beta,
            "ordinary_coefficient_bound": coefficient_bound,
            "input_digit_bound": max(maximum_a, maximum_b),
            "block_bits": 8 * width,
            "packed_a_bits": packed_a.bit_length(),
            "packed_b_bits": packed_b.bit_length(),
            "packed_product_bits": packed_product.bit_length(),
            "packed_product_bytes_allocated": raw_product_bytes,
            "squaring_same_object": a is b,
        }
    data = memoryview(packed_product.to_bytes(raw_product_bytes, "little"))
    del packed_a, packed_b, packed_product
    out = []
    for k in range(m):
        low = int.from_bytes(data[k * width:(k + 1) * width], "little")
        high = 0 if k == m - 1 else int.from_bytes(
            data[(k + m) * width:(k + m + 1) * width], "little")
        out.append(low + high - correction)
    if observations is not None:
        observation["output_max_integer_bits"] = max(abs(x).bit_length() for x in out)
        observation["output_nonzero_elements"] = sum(x != 0 for x in out)
        observations.append(observation)
    return out


def run(
    N: int,
    L: int,
    f: int,
    params: list[int],
    instrument: bool = False,
    method: str = "binary",
) -> dict:
    """Same successful-result timing/output contract as reference.run.

    Memory infeasibility returns status='memory_precheck_blocked' without
    constructing a field or giant array and without a completion time.
    Normal runs have no diagnostic counts/array scans.  Scans needed to choose
    an exact carry-free block width and offset correction are algorithm work
    and remain charged.  instrument=True timings are ineligible for aggregation.
    """
    if method not in ("binary", "kronecker", "kronecker_binary"):
        raise ValueError("This module implements Kronecker binary powering only")
    if N < 1 or L < 1:
        raise ValueError("N and L must be positive")
    if not params or any(not 0 < a < 1 << N for a in params):
        raise ValueError("At least one nonzero field parameter is required")
    start = perf_counter()
    resource_check = memory_precheck(N)
    if not resource_check["feasible_by_pointer_lower_bound"]:
        return {
            "status": "memory_precheck_blocked",
            "method": "convolution_kronecker_binary",
            "N": N, "L": L, "f": f, "params": params,
            "outputs": None, "compute_seconds": None, "cold_seconds": None,
            "preprocess_seconds": None, "query_seconds": None,
            "memory_precheck": resource_check,
            "timing_eligible": False,
            "note": "No full-sized arrays allocated; field validity not evaluated after resource rejection.",
        }
    own_trace = instrument and not tracemalloc.is_tracing()
    if own_trace:
        tracemalloc.start()
    observations = [] if instrument else None
    if not reference.irreducible(N, f):
        raise ValueError("The supplied modulus is not irreducible of degree N")
    primitive, exponents, logarithms = reference._cyclic_index(N, f)
    after_index = perf_counter()
    base = [1 - 2 * reference.absolute_trace(a, N, f) for a in exponents]
    del exponents
    after_character = perf_counter()
    m = (1 << N) - 1
    table = [1] + [0] * (m - 1)
    exponent = L
    while exponent:
        if exponent & 1:
            table = cyclic_kronecker(table, base, observations)
        exponent >>= 1
        if exponent:
            base = cyclic_kronecker(base, base, observations)
    del base
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
        "status": "ok",
        "method": "convolution_kronecker_binary",
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
        "memory_precheck": resource_check,
        "table_is_all_parameters": True,
        "reusable_objects": ["logarithm_index", "full_convolution_table"],
        "integer_engine": "CPython built-in arbitrary-precision int multiplication",
    }
    if instrument:
        current_bytes, peak_bytes = tracemalloc.get_traced_memory()
        if own_trace:
            tracemalloc.stop()
        result["instrumentation"] = {
            "tracemalloc_current_before_census_bytes": current_bytes,
            "tracemalloc_peak_before_census_bytes": peak_bytes,
            "tracemalloc_scope": "this_run" if own_trace else "external_tracer",
            "retained_array_slots": len(logarithms) + len(table),
            "retained_arrays_python_bytes_deduplicated": reference._retained_python_bytes((logarithms, table)),
            "retained_arrays": {
                "logarithm_index": reference._array_summary(logarithms),
                "full_convolution_table": reference._array_summary(table),
            },
            "packed_integer_products": len(observations),
            "packing_records": observations,
            "maximum_observed_value_bits": max(r["output_max_integer_bits"] for r in observations),
            "note": "This traced run includes observation metadata. For allocation-only peaks, externally trace run(instrument=False). Census scratch is excluded from the peak sampled above.",
        }
    return result


def check() -> dict:
    """Small, untimed-as-performance validation of the new exact arithmetic."""
    import hashlib
    import random

    rng = random.Random(20260930)
    cases = [([0], [0]), ([-2] * 3, [-5] * 3), ([1, 0, 0], [-1, 7, -9]),
             ([0] * 5, [2, -3, 4, -5, 6]), ([0, 0], [1 << 100, -(1 << 99)])]
    for length in (1, 2, 5, 15):
        for bits in (2, 8, 61):
            cases.append(([rng.randrange(-(1 << bits), 1 << bits) for _ in range(length)],
                          [rng.randrange(-(1 << bits), 1 << bits) for _ in range(length)]))
    for a, b in cases:
        assert cyclic_kronecker(a, b) == reference.cyclic_general(a, b)
        assert cyclic_kronecker(a, a) == reference.cyclic_general(a, a)
    # The frozen literal checker imports 'convolution'.  Restrict that import
    # to its own frozen reference module while executing this validation only.
    prior_convolution = sys.modules.get("convolution")
    sys.modules["convolution"] = reference
    try:
        literal = _load_file("_qcb_literal_reference", REFERENCE_CODE / "check_convolution.py")
    finally:
        if prior_convolution is None:
            del sys.modules["convolution"]
        else:
            sys.modules["convolution"] = prior_convolution
    checked = []
    for N, f, ranks in ((2, 0x7, (1, 2, 3, 4)), (3, 0xB, (2, 3)), (4, 0x13, (2, 3))):
        params = list(range(1, 1 << N))
        for L in ranks:
            expected = literal.literal_values(N, L, f, params)
            result = run(N, L, f, params)
            assert result["outputs"] == expected
            checked.append({"N": N, "L": L, "params": params, "outputs": expected,
                            "all_parameters_equal_literal": True})
    instrumented = run(4, 4, 0x13, [2, 1, 3, 5], instrument=True)
    assert instrumented["outputs"] == reference.run(4, 4, 0x13, [2, 1, 3, 5])["outputs"]
    assert instrumented["instrumentation"]["packed_integer_products"] == 3
    blocked = [run(N, 4, 0, [2]) for N in (32, 64)]
    assert all(row["status"] == "memory_precheck_blocked" and row["compute_seconds"] is None
               for row in blocked)
    return {
        "status": "passed",
        "purpose": "bounded arithmetic validation; no performance timings retained",
        "implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "signed_and_square_cases": len(cases),
        "literal_cases": checked,
        "instrumentation_count_checked": True,
        "memory_precheck_cases": blocked,
        "scope": "Exactness argument plus these finite cases; no empirical large-field claim.",
    }


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=4)
    parser.add_argument("--L", type=int, default=4)
    parser.add_argument("--f", type=lambda s: int(s, 0))
    parser.add_argument("--params", type=int, nargs="+", default=[2, 1, 3, 5])
    parser.add_argument("--instrument", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not __debug__ and args.check:
        parser.error("Do not disable assertions in --check")
    if args.check:
        record = check()
    else:
        modulus = args.f if args.f is not None else FIELD_POLYNOMIALS.get(args.N, 0)
        record = run(args.N, args.L, modulus, args.params, args.instrument)
    serialized = json.dumps(record, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(serialized, end="", flush=True)
