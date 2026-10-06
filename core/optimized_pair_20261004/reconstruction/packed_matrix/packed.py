"""D6 exact Gaussian quotient-ring matrix products via public fmpz_mat.

No scientific operation is performed on import. See PROTOCOL.md and REVIEW.md.
The helpers use Python integers/bytes; native matrices never expose private
object layouts. All packing is exact, with no loss of precision.
"""
from struct import calcsize
import sys
from time import perf_counter

from flint import fmpz_mat, fmpz_mod_poly, fmpz_mod_poly_ctx


IMPLEMENTATION = "packed_gaussian_matrix_d6_v1"
_LIMIT = min(sys.maxsize, (1 << (8 * calcsize("l") - 1)) - 1)
_POINTER_BYTES = calcsize("P")
_STAGES = ("check", "pack", "native", "decode", "reduce")


def _positive(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(name + " must be a positive Python integer")
    if value > _LIMIT:
        raise OverflowError(name + " exceeds the Python/FLINT size limit")
    return value


def _size_product(left, right, name):
    if left > _LIMIT // right:
        raise OverflowError(name + " exceeds the Python/FLINT size limit")
    return left * right


def _buffer_layout(length, field_bytes):
    _positive(length, "length")
    _positive(field_bytes, "field_bytes")
    field_bits = _size_product(field_bytes, 8, "field bit width")
    byte_length = _size_product(length, field_bytes, "byte buffer length")
    # The packed bit length also enters FLINT's signed length arithmetic.
    _size_product(length, field_bits, "packed bit length")
    return field_bits, byte_length


def pack_coefficients(coefficients, P, d, field_bytes):
    """Encode degree<d canonical Python integer coefficients into one int.

    Omitted high zero coefficients are permitted; oversized lists are not.
    The complete capacity is checked before building any byte buffer.
    """
    _positive(P, "P")
    field_bits, _ = _buffer_layout(d, field_bytes)
    if field_bits < P:
        raise ValueError("field cannot hold a canonical P-bit coefficient")
    if not isinstance(coefficients, (list, tuple)) or len(coefficients) > d:
        raise ValueError("expected at most d coefficients")
    for value in coefficients:
        if type(value) is not int or value < 0 or value.bit_length() > P:
            raise ValueError("noncanonical input coefficient")
    # Missing high fields are zero. One join avoids repeatedly copying a
    # growing packed integer; from_bytes is independent of decimal limits.
    data = b"".join(value.to_bytes(field_bytes, "little") for value in coefficients)
    return int.from_bytes(data, "little")


def unpack_coefficients(value, length, field_bytes, signed=False):
    """Recover a fixed number of unsigned or strictly balanced coefficients.

    signed=True uses the whole integer's sign, plus a one-bit inter-field
    carry. It does not interpret ordinary two's-complement fields separately.
    A half-base coefficient and a nonzero final carry are rejected.
    """
    if type(value) is not int or type(signed) is not bool:
        raise ValueError("expected a Python integer and boolean signed flag")
    field_bits, byte_length = _buffer_layout(length, field_bytes)
    if not signed and value < 0:
        raise ArithmeticError("negative integer in unsigned unpack")
    magnitude = abs(value)
    if magnitude.bit_length() > length * field_bits:
        raise ArithmeticError("packed integer exceeds the declared length")
    data = magnitude.to_bytes(byte_length, "little")
    view = memoryview(data)
    if not signed:
        return [int.from_bytes(view[j:j + field_bytes], "little")
                for j in range(0, byte_length, field_bytes)]
    base = 1 << field_bits
    half = base >> 1
    sign = -1 if value < 0 else 1
    carry = 0
    out = []
    for j in range(0, byte_length, field_bytes):
        digit = int.from_bytes(view[j:j + field_bytes], "little") + carry
        if digit >= half:
            coefficient, carry = digit - base, 1
        else:
            coefficient, carry = digit, 0
        if not -half < coefficient < half:
            raise ArithmeticError("half-base coefficient in balanced unpack")
        out.append(sign * coefficient)
    if carry:
        raise ArithmeticError("nonzero final carry in balanced unpack")
    return out


class PackedMultiplier:
    """One query's finite ring and paid multiplication statistics.

    calls counts attempts; successful_calls + failed_calls == calls.
    Branch counters count successful returns. native_products/mod_g_calls
    count completed primitive calls, including work before a later failure.
    All time totals include failed attempts; constructor checks are charged
    to check and total. Native/input temporaries are released before total
    stops; the returned matrix's lifetime belongs to the enclosing query.
    """

    def __init__(self, R):
        started = perf_counter()
        self._seconds = {stage: 0.0 for stage in _STAGES}
        self._seconds["total"] = 0.0
        self._counts = dict(calls=0, successful_calls=0, failed_calls=0,
                            zero_calls=0, real_calls=0, one_imag_calls=0,
                            gaussian_calls=0, native_products=0, mod_g_calls=0,
                            input_coefficients=0, packed_entries=0,
                            packed_coefficients=0, decoded_coefficients=0,
                            max_field_bits=0)
        self.P = _positive(getattr(R, "e", None), "ring precision")
        self.d = _positive(getattr(R, "d", None), "ring degree")
        modulus = getattr(R, "mod", None)
        if (type(modulus) is not int or modulus <= 1
                or modulus.bit_length() != self.P + 1
                or modulus & (modulus - 1)):
            raise ValueError("ring modulus must equal 2**P")
        ctx, g = getattr(R, "ctx", None), getattr(R, "g", None)
        if not isinstance(ctx, fmpz_mod_poly_ctx) or not isinstance(g, fmpz_mod_poly):
            raise ValueError("expected public modular-polynomial ring objects")
        if int(ctx.modulus()) != modulus or int(g.modulus()) != modulus:
            raise ValueError("ring context modulus mismatch")
        if g.degree() != self.d or int(g[self.d]) != 1:
            raise ValueError("expected a monic degree-d modulus polynomial")
        self.mod, self.mask, self.ctx, self.g = modulus, modulus - 1, ctx, g
        self.zero = (ctx([]), ctx([]))
        self._layout(1, 1, 1)
        elapsed = perf_counter() - started
        self._seconds["check"] = elapsed
        self._seconds["total"] = elapsed

    def statistics_snapshot(self):
        return dict(implementation=IMPLEMENTATION, bits=self.P, degree=self.d,
                    **self._counts, seconds=dict(self._seconds))

    def _layout(self, m, k, n):
        for value, name in ((m, "rows"), (k, "inner dimension"), (n, "columns")):
            _positive(value, name)
        for left, right in ((m, k), (k, n), (m, n)):
            entries = _size_product(left, right, "matrix entry count")
            _size_product(entries, _POINTER_BYTES, "matrix entry storage")
        kd = _size_product(k, self.d, "inner convolution size")
        guard = (kd - 1).bit_length() + 1
        if self.P > (_LIMIT - guard - 7) // 2:
            raise OverflowError("packing width exceeds the size limit")
        field_bytes = (2 * self.P + guard + 7) // 8
        if self.d > (_LIMIT + 1) // 2:
            raise OverflowError("convolution length exceeds the size limit")
        length = 2 * self.d - 1
        field_bits, _ = _buffer_layout(length, field_bytes)
        _size_product(m * n, length, "decoded coefficient count")
        return field_bytes, field_bits, length

    @staticmethod
    def _shape(matrix):
        if not isinstance(matrix, (list, tuple)) or not matrix:
            raise ValueError("expected a nonempty matrix")
        if not isinstance(matrix[0], (list, tuple)) or not matrix[0]:
            raise ValueError("expected nonempty matrix rows")
        columns = len(matrix[0])
        if any(not isinstance(row, (list, tuple)) or len(row) != columns
               for row in matrix):
            raise ValueError("ragged matrix")
        return len(matrix), columns

    def _coefficients(self, matrix):
        real, imag = [], []
        nonzero, has_imag = False, False
        for row in matrix:
            for value in row:
                if not isinstance(value, (list, tuple)) or len(value) != 2:
                    raise ValueError("expected real/imaginary polynomial pairs")
                for part, target in enumerate((real, imag)):
                    polynomial = value[part]
                    if not isinstance(polynomial, fmpz_mod_poly):
                        raise ValueError("matrix component is not a modular polynomial")
                    if int(polynomial.modulus()) != self.mod:
                        raise ValueError("matrix component context mismatch")
                    if polynomial.degree() >= self.d:
                        raise ValueError("matrix component degree is not below d")
                    coefficients = [int(polynomial[j]) for j in range(len(polynomial))]
                    self._counts["input_coefficients"] += len(coefficients)
                    if any(c < 0 or c >= self.mod for c in coefficients):
                        raise ValueError("noncanonical matrix coefficient")
                    active = any(coefficients)
                    nonzero |= active
                    if part:
                        has_imag |= active
                    target.append(coefficients)
        return real, imag, nonzero, has_imag

    def _pack_matrix(self, coefficients, rows, columns, width):
        values = []
        for polynomial in coefficients:
            values.append(pack_coefficients(polynomial, self.P, self.d, width))
            self._counts["packed_entries"] += 1
            self._counts["packed_coefficients"] += len(polynomial)
        return fmpz_mat(rows, columns, values)

    def _product(self, left, right):
        product = left * right
        self._counts["native_products"] += 1
        return product

    def _decode_matrix(self, matrix, m, n, length, width, signed):
        out = []
        for i in range(m):
            for j in range(n):
                coefficients = unpack_coefficients(int(matrix[i, j]), length, width, signed)
                self._counts["decoded_coefficients"] += length
                out.append([c & self.mask for c in coefficients])
        return out

    def multiply(self, A, B):
        started = phase_started = perf_counter()
        phase = "check"
        self._counts["calls"] += 1
        work = {}
        left = right = real = imag = None
        branch = None
        success = False

        def next_phase(name):
            nonlocal phase, phase_started
            now = perf_counter()
            self._seconds[phase] += now - phase_started
            phase, phase_started = name, now

        try:
            m, k = self._shape(A)
            inner, n = self._shape(B)
            if inner != k:
                raise ValueError("incompatible matrix dimensions")
            width, field_bits, length = self._layout(m, k, n)
            left, right = self._coefficients(A), self._coefficients(B)
            self._counts["max_field_bits"] = max(self._counts["max_field_bits"], field_bits)
            if not left[2] or not right[2]:
                branch = "zero_calls"
                next_phase("reduce")
                result = [[self.zero for _ in range(n)] for _ in range(m)]
            else:
                has_ai, has_bi = left[3], right[3]
                branch = ("gaussian_calls" if has_ai and has_bi else
                          "one_imag_calls" if has_ai or has_bi else "real_calls")
                next_phase("pack")
                work["Ar"] = self._pack_matrix(left[0], m, k, width)
                work["Br"] = self._pack_matrix(right[0], k, n, width)
                if has_ai:
                    work["Ai"] = self._pack_matrix(left[1], m, k, width)
                if has_bi:
                    work["Bi"] = self._pack_matrix(right[1], k, n, width)
                next_phase("native")
                work["U"] = self._product(work["Ar"], work["Br"])
                if has_ai and has_bi:
                    work["V"] = self._product(work["Ai"], work["Bi"])
                    work["sumA"] = work["Ar"] + work["Ai"]
                    work["sumB"] = work["Br"] + work["Bi"]
                    work["W"] = self._product(work["sumA"], work["sumB"])
                    work["E"] = work["U"] - work["V"]
                    # W may contain carries. Keep its complete exact value
                    # until both subtractions have formed the cross term.
                    work["I"] = work["W"] - work["U"] - work["V"]
                else:
                    work["E"] = work["U"]
                    if has_ai:
                        work["I"] = self._product(work["Ai"], work["Br"])
                    elif has_bi:
                        work["I"] = self._product(work["Ar"], work["Bi"])
                next_phase("decode")
                real = self._decode_matrix(work["E"], m, n, length, width, has_ai and has_bi)
                if has_ai or has_bi:
                    imag = self._decode_matrix(work["I"], m, n, length, width, False)
                next_phase("reduce")
                result = []
                for i in range(m):
                    row = []
                    for j in range(n):
                        position = i * n + j
                        re = self.ctx(real[position]) % self.g
                        self._counts["mod_g_calls"] += 1
                        if imag is None:
                            im = self.zero[1]
                        else:
                            im = self.ctx(imag[position]) % self.g
                            self._counts["mod_g_calls"] += 1
                        row.append((re, im))
                    result.append(row)
            success = True
            return result
        finally:
            self._seconds[phase] += perf_counter() - phase_started
            # Release native matrices and coefficient work before total stops.
            # Output matrix ownership transfers to the query on success.
            work.clear()
            left = right = real = imag = None
            if success:
                self._counts["successful_calls"] += 1
                self._counts[branch] += 1
            else:
                self._counts["failed_calls"] += 1
            self._seconds["total"] += perf_counter() - started
