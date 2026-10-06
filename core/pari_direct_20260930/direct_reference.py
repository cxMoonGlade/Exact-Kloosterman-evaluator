"""Independent direct binary character sums and polynomial irreducibility."""

MODULI = {2: 7, 4: 19, 6: 67, 8: 283, 12: 4179}


def remainder(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    return a


def multiply(a, b, modulus):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
    return remainder(result, modulus)


def power(a, exponent, modulus):
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, a, modulus)
        a = multiply(a, a, modulus)
        exponent >>= 1
    return result


def irreducible(N, modulus):
    if modulus.bit_length() != N + 1 or modulus % 2 == 0:
        return False
    frobenius = 2
    for degree in range(1, N + 1):
        frobenius = multiply(frobenius, frobenius, modulus)
        if degree <= N // 2 and gcd(frobenius ^ 2, modulus) != 1:
            return False
    return frobenius == 2


def character(a, N, modulus):
    trace = 0
    for _ in range(N):
        trace ^= a
        a = multiply(a, a, modulus)
    assert trace in (0, 1)
    return 1 - 2 * trace


def sums(N, parameters):
    modulus = MODULI[N]
    q = 1 << N
    inverses = [0] + [power(x, q - 2, modulus) for x in range(1, q)]
    characters = [character(x, N, modulus) for x in range(q)]
    result = {}
    for a in parameters:
        rank2 = sum(characters[x ^ multiply(a, inverses[x], modulus)]
                    for x in range(1, q))
        rank3 = sum(characters[x ^ y ^ multiply(
            a, multiply(inverses[x], inverses[y], modulus), modulus)]
            for x in range(1, q) for y in range(1, q))
        result[(2, a)] = rank2
        result[(3, a)] = rank3
    return result
