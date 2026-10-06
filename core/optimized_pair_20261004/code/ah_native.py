"""Finite AH acquisition through two native FLINT multiple-RHS solves.

See NATIVE_API_AH.md and NATIVE_API_REVIEW.md for composite-modulus LU
applicability and the finite representative chosen for each column of V.
This dense branch does not claim the structured acquisition complexity.
No scientific computation occurs on import; all contexts are local/cold.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
from time import perf_counter

import flint
from flint import fmpz_mod_ctx, fmpz_mod_mat, fmpz_mod_poly_ctx


SCALARS = (Path(__file__).resolve().parents[2] /
           "specialization_optimization_20260930/reference/code/modular_scalars.py")


def _scalar_data(bits):
    spec = importlib.util.spec_from_file_location("native_pair_frozen_c0", SCALARS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.scalar_data(bits)


def _check_native_class(rows, bits):
    """Restrict the native first-nonzero pivot chain to odd pivots."""
    if type(bits) is not int or bits < 1:
        raise ValueError("positive integer modulus precision required")
    n = len(rows)
    if not n or any(len(row) != n for row in rows):
        raise ValueError("nonempty square matrix required")
    if any(not (rows[i][i] & 1) or
           any(rows[i][j] & 1 for j in range(i)) for i in range(n)):
        raise ValueError("native LU requires unit-upper-triangular mod-2 input")
    if flint.__version__ != "0.9.0" or flint.__FLINT_VERSION__ != "3.6.0":
        raise RuntimeError("unreviewed native LU version")


def _residual_valid(rows, columns, rhs, bits):
    """Independent literal modular residual, diagnostic only."""
    mask = (1 << bits) - 1
    return all((sum(a*x for a, x in zip(row, col)) - expected[i]) & mask == 0
               for col, expected in zip(columns, rhs)
               for i, row in enumerate(rows))


def _solve(rows, rhs_columns, bits, *, audit=False):
    _check_native_class(rows, bits)
    n = len(rows)
    if not rhs_columns or any(len(col) != n for col in rhs_columns):
        raise ValueError("nonempty RHS columns must match the system dimension")
    # Both operands are constructed here with the very same context. The public
    # wrapper does not itself certify the contexts of supplied operands agree.
    ctx = fmpz_mod_ctx(1 << bits)
    matrix = fmpz_mod_mat(rows, ctx)
    rhs = fmpz_mod_mat([list(row) for row in zip(*rhs_columns)], ctx)
    solution = matrix.solve(rhs)
    columns = [[int(solution[i, j]) for i in range(n)]
               for j in range(len(rhs_columns))]
    if audit and not _residual_valid(rows, columns, rhs_columns, bits):
        raise ArithmeticError("native multiple-RHS solve residual failed")
    return columns


def acquire(N, L, P=None, *, audit=False):
    if type(N) is not int or type(L) is not int or N < 1 or L < 2:
        raise ValueError("positive N and L>=2 required")
    P = N*(L-1)+2 if P is None else P
    if type(P) is not int or P < N*(L-1)+2:
        raise ValueError("output precision below exact recovery requirement")
    Q = P+L-1
    mod, mask = 1 << Q, (1 << Q)-1
    stages, counts = {}, {}
    started = perf_counter()

    t = perf_counter()
    gamma, z, phi_values, b = _scalar_data(Q)
    if len(phi_values) != Q or len(b) < 2*Q:
        raise ValueError("incompatible cold scalar-constructor output")
    stages["scalar_construction"] = perf_counter()-t

    t = perf_counter()
    poly_ctx = fmpz_mod_poly_ctx(mod)
    unit = poly_ctx([1])
    phi = poly_ctx(phi_values)
    del phi_values, z
    products = 0

    def multiply(a, b):
        nonlocal products
        products += 1
        return a.mul_low(b, Q)

    def powers(base, count):
        result = [unit]
        if count:
            result.append(base)
            for _ in range(2, count+1):
                result.append(multiply(result[-1], base))
        return result

    U = [[0]*Q for _ in range(Q)]
    u = unit
    for k in range(Q):
        column = [int(u[r]) for r in range(Q)]
        for r, value in enumerate(column):
            U[r][k] = value
        if k+1 < Q:
            product = multiply(u, phi)
            u = poly_ctx([(r*x+int(product[r])) & mask
                          for r, x in enumerate(column)])
    counts["first_transform_products"] = products
    del u, column, phi
    stages["first_transform"] = perf_counter()-t

    t = perf_counter()
    gamma_powers = [1]
    for _ in range(Q):
        gamma_powers.append((gamma_powers[-1]*gamma) & mask)
    required = (Q+L-2)//L+1
    rhs_u = [[(1 << k) if r == k else 0 for r in range(Q)]
             for k in range(required)]
    for offset in (-1, 0, 1):
        rhs_u.append([0 if 2*k+offset < 0 else
                      (gamma_powers[k+(offset == 1)]*b[2*k+offset]) & mask
                      for k in range(Q)])
    del b
    stages["first_rhs_construction"] = perf_counter()-t

    t = perf_counter()
    solved = _solve(U, rhs_u, Q, audit=audit)
    hats = [poly_ctx(col) for col in solved[:required]]
    hm, h0, hp = (poly_ctx(col) for col in solved[required:])
    del U, rhs_u, solved
    counts.update(native_factorizations=1, first_solve_rhs=required+3,
                  first_input_matrix_slots=Q*Q, first_rhs_slots=Q*(required+3))
    stages["first_batch_solve"] = perf_counter()-t

    t = perf_counter()
    V = [[0]*Q for _ in range(Q)]
    before = products

    def write_normalized(poly, k):
        divisor_mask = (1 << k)-1
        for r in range(Q):
            value = int(poly[r])
            if value & divisor_mask:
                raise ArithmeticError("finite V column division is not exact")
            # This is the proved finite lift: Q-k known low bits, high bits 0.
            V[r][k] = value >> k

    running = unit
    for j in range(min(L, Q)):
        if j == 1:
            running = hats[1]
        elif j > 1:
            running = multiply(running, hats[1])
        write_normalized(running, j)
    for d in range(1, (Q-1)//L+1):
        apowers = powers(hats[d], L)
        running = unit
        for j in range(min(L, Q-L*d)):
            if j == 1:
                running = hats[d+1]
            elif j > 1:
                running = multiply(running, hats[d+1])
            weighted = apowers[L-j] if j == 0 else multiply(apowers[L-j], running)
            write_normalized(weighted, L*d+j)
        del apowers, weighted
    del hats, running
    counts["second_transform_products"] = products-before
    counts["exact_column_shifts"] = Q*Q
    stages["second_transform"] = perf_counter()-t

    t = perf_counter()
    before = products
    h0powers, hppowers = powers(h0, L), powers(hp, L)
    hm_power = unit
    rhs_v = []
    for i in range(L):
        if i == 1:
            hm_power = hm
        elif i > 1:
            hm_power = multiply(hm_power, hm)
        s0 = h0powers[L-i] if i == 0 else multiply(h0powers[L-i], hm_power)
        s1 = hppowers[L] if i == 0 else multiply(hppowers[L-i], h0powers[i])
        rhs_v.append([int(s0[k]) for k in range(Q)])
        rhs_v.append([(int(s1[k])*gamma_powers[i]) & mask for k in range(Q)])
    counts["second_rhs_products"] = products-before
    counts["polynomial_products"] = products
    del hm, h0, hp, h0powers, hppowers, hm_power, s0, s1
    stages["second_rhs_construction"] = perf_counter()-t

    t = perf_counter()
    solved = _solve(V, rhs_v, Q, audit=audit)
    del V, rhs_v
    counts.update(native_factorizations=2, second_solve_rhs=2*L,
                  second_input_matrix_slots=Q*Q, second_rhs_slots=Q*2*L)
    stages["second_batch_solve"] = perf_counter()-t

    t = perf_counter()
    out = [[[0]*((Q-1-j)//L+1) for _ in range(L)] for j in range(L)]
    for i in range(L):
        for h in range(2):
            for k, value in enumerate(solved[2*i+h]):
                if value & ((1 << min(Q, k+L*h))-1):
                    raise ArithmeticError("native solution escaped the legal filtered module")
                d, j = divmod(k, L)
                degree = 2*d+h
                coefficient = (value*gamma_powers[L*d]) & mask
                if L*degree+j < Q:
                    out[j][i][degree] = coefficient
                elif coefficient:
                    raise ArithmeticError("nonzero coefficient outside proved support")
    stages["supported_output"] = perf_counter()-t
    elapsed = perf_counter()-started
    return dict(coefficients=out, gamma=gamma, P=P, Q=Q, stages=stages,
                counts=counts, acquisition_seconds=elapsed,
                stage_gap_seconds=elapsed-sum(stages.values()),
                precision=dict(output=P, acquisition=Q, native_modulus_bits=Q,
                               column_lift="standard representative after exact division"),
                diagnostics=dict(enabled=bool(audit),
                                 full_system_residuals=2 if audit else 0),
                native=dict(python_flint=flint.__version__, flint=flint.__FLINT_VERSION__),
                layout="coefficients[target_row][source_column][degree]",
                algorithm="native_flint_two_batch_solves")


def validate_components():
    """Small class, residual, and frozen-finite-array checks; no timing claim."""
    from ah_regression import validate_regression
    regression = validate_regression(lambda N, L: acquire(N, L, audit=True))
    rows = [[3, 7, 2], [2, 5, 4], [6, 2, 7]]
    rhs = [[1, 2, 3], [5, 7, 11], [0, 0, 0]]
    solution = _solve(rows, rhs, 9, audit=True)
    changed = [col.copy() for col in solution]
    changed[1][0] ^= 1
    if _residual_valid(rows, changed, rhs, 9):
        raise AssertionError("changed solution escaped the independent residual check")
    rejected = 0
    for bad in ([[2, 0], [0, 1]], [[1, 0], [1, 1]], [[1, 2, 3], [0, 1]]):
        try:
            _check_native_class(bad, 9)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("unsupported native LU input class was accepted")
    return dict(status="PASS", regression=regression, residual_mutation="DETECTED",
                invalid_matrix_classes_rejected=rejected,
                native_solve_residuals="both complete U and V systems in all regression cases")
