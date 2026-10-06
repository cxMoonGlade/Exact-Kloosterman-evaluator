"""Shared exact parameter presentation and balanced integral flag.

Premises and resource scope: ../PROTOCOL.md and PROTOCOL_REVIEW.md.
All polynomial products/remainders use FLINT. Matrices are small classical
blocks; no theorem-level fast matrix multiplication claim is made here.
There is no scientific computation on import.
"""
from time import perf_counter
from flint import fmpz_poly, fmpz_mod_poly_ctx


def bits(a):
    return [(a >> j) & 1 for j in range(max(1, a.bit_length()))]


def gf_reduce(a, f):
    n = f.bit_length() - 1
    while a.bit_length() > n:
        a ^= f << (a.bit_length() - n - 1)
    return a


def gf_mul(a, b, f):
    out = 0
    while b:
        if b & 1:
            out ^= a
        a <<= 1
        b >>= 1
    return gf_reduce(out, f)


def gf_square(a, f):
    out = sum(((a >> j) & 1) << (2*j) for j in range(a.bit_length()))
    return gf_reduce(out, f)


def gf_gcd(a, b):
    while b:
        a, b = b, gf_reduce(a, b)
    return a


def validate_input(N, L, f, parameters):
    if N < 1 or L < 2 or f.bit_length() != N+1:
        raise ValueError("expected N>=1, L>=2, monic degree-N binary f")
    # Rabin irreducibility test, including both valid degree-one moduli.
    x = gf_reduce(2, f)
    y = x
    for j in range(1, N+1):
        y = gf_square(y, f)
        if j <= N//2 and gf_gcd(y ^ x, f) != 1:
            raise ValueError("reducible field modulus")
    if y != x or any(a <= 0 or a >= 1 << N for a in parameters):
        raise ValueError("invalid field modulus or nonzero parameter")


def parameter_polynomial(a, f):
    """Actual Frobenius orbit; native packed bivariate product tree over F2."""
    orbit, value = [], a
    while value not in orbit:
        orbit.append(value)
        value = gf_square(value, f)
    assert value == a
    N, d = f.bit_length()-1, len(orbit)
    assert N % d == 0
    ctx = fmpz_mod_poly_ctx(2)
    width = 2*N

    def multiply(p, q):
        def pack(v):
            coefficients = [0]*(width*(len(v)-1)+N)
            for j, c in enumerate(v):
                for i in range(N):
                    coefficients[width*j+i] = (c >> i) & 1
            return ctx(coefficients)
        pq = pack(p)*pack(q)
        return [gf_reduce(sum(int(pq[width*j+i]) << i
                              for i in range(width-1)), f)
                for j in range(len(p)+len(q)-1)]

    level = [[root, 1] for root in orbit]
    while len(level) > 1:
        level = [multiply(level[j], level[j+1]) if j+1 < len(level)
                 else level[j] for j in range(0, len(level), 2)]
    assert all(c in (0, 1) for c in level[0])
    return level[0], orbit


def contraction_vector(b, e, operator):
    """Solve x-T(x)=b; T raises binary valuation, at actual local widths."""
    mask = (1 << e)-1
    b = [c & mask for c in b]
    if e == 1:
        return b
    m = (e+1)//2
    lo = contraction_vector(b, m, operator)
    tlo = operator(lo, e)
    residual = [(c-x+t) & mask for c, x, t in zip(b, lo, tlo)]
    assert all(c % (1 << m) == 0 for c in residual)
    hi = contraction_vector([c >> m for c in residual], e-m, operator)
    return [x+(y << m) for x, y in zip(lo, hi)]


def graeffe(g, e):
    sign = (-1)**(len(g)-1)
    p = fmpz_poly(g)*fmpz_poly([c*(-1)**j for j, c in enumerate(g)])
    mask = (1 << e)-1
    return [(sign*int(p[2*j])) & mask for j in range(len(g))]


def teich_modulus(fa, e):
    g, known, d = list(fa), 1, len(fa)-1
    while known < e:
        target = min(2*known, e)
        mask = (1 << target)-1
        gg = graeffe(g, target)
        rhs = [(gg[j]-g[j]) & mask for j in range(d)]

        def operator(h, width):
            prefix = (1 << width)-1
            p = fmpz_poly([c & prefix for c in h])
            q = fmpz_poly([(c & prefix)*(-1)**j for j, c in enumerate(g)])
            product = p*q
            return [(2*(-1)**d*int(product[2*j])) & prefix for j in range(d)]

        h = contraction_vector(rhs, target, operator)
        g = [(g[j]+h[j]) & mask for j in range(d)]+[1]
        assert graeffe(g, target) == g
        known = target
    return g


class Ring:
    """S_e[i], with i fixed by Frobenius, represented by two native polynomials."""
    def __init__(self, presentation, e):
        self.parent, self.e = presentation, e
        self.mod = 1 << e
        self.mask = self.mod-1
        self.ctx = fmpz_mod_poly_ctx(self.mod)
        self.g = self.ctx(presentation.g)
        self.d = len(presentation.g)-1
        self.zero = (self.ctx([]), self.ctx([]))
        self.one = (self.ctx([1]), self.ctx([]))
        self.xi = None

    def elem(self, real=(), imag=()):
        return self.ctx(list(real)) % self.g, self.ctx(list(imag)) % self.g

    def reduce(self, a):
        return self.elem([int(c) for c in a[0]], [int(c) for c in a[1]])

    def add(self, a, b):
        return a[0]+b[0], a[1]+b[1]

    def sub(self, a, b):
        return a[0]-b[0], a[1]-b[1]

    def scale(self, a, c):
        return a[0]*c, a[1]*c

    def mul(self, a, b):
        if not a[1] and not b[1]:
            return (a[0]*b[0]) % self.g, self.ctx([])
        return ((a[0]*b[0]-a[1]*b[1]) % self.g,
                (a[0]*b[1]+a[1]*b[0]) % self.g)

    def base_inverse(self, a):
        low = self.parent.ring(1)
        y = low.ctx([int(c) for c in a]).inverse_mod(low.g)
        width = 1
        while width < self.e:
            width = min(2*width, self.e)
            rr = self.parent.ring(width)
            yy = rr.ctx([int(c) for c in y])
            aa = rr.ctx([int(c) for c in a])
            y = (yy*(2-aa*yy)) % rr.g
        return self.ctx([int(c) for c in y])

    def inv(self, a):
        denominator = (a[0]*a[0]+a[1]*a[1]) % self.g
        z = self.base_inverse(denominator)
        return (a[0]*z) % self.g, (-a[1]*z) % self.g

    def pow(self, a, n):
        out = self.one
        while n:
            if n & 1:
                out = self.mul(out, a)
            n >>= 1
            if n:
                a = self.mul(a, a)
        return out

    def sigma(self, a):
        if self.d == 1:
            return a
        def square_argument(p):
            coeff = [0]*(max(0, 2*len(p)-1))
            for j, c in enumerate(p):
                coeff[2*j] = int(c)
            return self.ctx(coeff) % self.g
        return square_argument(a[0]), square_argument(a[1])

    def tau(self, a):
        if self.d == 1:
            return a
        if self.xi is None:
            if self.parent.xi is None:
                self.parent.prepare_tau()
            self.xi = self.ctx(self.parent.xi) % self.g
        def inverse_argument(p):
            even = self.ctx([int(p[j]) for j in range(0, len(p), 2)])
            odd = self.ctx([int(p[j]) for j in range(1, len(p), 2)])
            return (even+self.xi*odd) % self.g
        return inverse_argument(a[0]), inverse_argument(a[1])

    def shift_down(self, a, s):
        assert all(int(c) % (1 << s) == 0 for p in a for c in p)
        return self.elem([int(c) >> s for c in a[0]],
                         [int(c) >> s for c in a[1]])


class Presentation:
    def __init__(self, fa, e):
        self.g = teich_modulus(fa, e)
        self.e, self.cache, self.xi = e, {}, None
        self.prepare_tau()

    def ring(self, e):
        assert 1 <= e <= self.e
        if e not in self.cache:
            self.cache[e] = Ring(self, e)
        return self.cache[e]

    def prepare_tau(self):
        if self.xi is not None:
            return
        R = self.ring(self.e)
        g0 = R.ctx(self.g[::2])
        g1 = R.ctx(self.g[1::2])
        self.xi = [int(c) for c in (-g0*R.base_inverse(g1)) % R.g]


def mzero(R, n, m):
    return [[R.zero for _ in range(m)] for _ in range(n)]


def identity(R, n):
    return [[R.one if i == j else R.zero for j in range(n)] for i in range(n)]


def mreduce(R, A):
    return [[R.reduce(v) for v in row] for row in A]


def madd(R, A, B):
    return [[R.add(x, y) for x, y in zip(a, b)] for a, b in zip(A, B)]


def msub(R, A, B):
    return [[R.sub(x, y) for x, y in zip(a, b)] for a, b in zip(A, B)]


def mscale(R, A, c):
    return [[R.scale(x, c) for x in row] for row in A]


def mmul(R, A, B):
    n, k, m = len(A), len(B), len(B[0])
    C = mzero(R, n, m)
    for i in range(n):
        for j in range(m):
            for t in range(k):
                C[i][j] = R.add(C[i][j], R.mul(A[i][t], B[t][j]))
    return C


def minverse(R, A):
    n = len(A)
    low = R.parent.ring(1)
    aa = mreduce(low, A)
    # Every required block is lower triangular with diagonal one mod 2.
    assert all(aa[i][j] == (low.one if i == j else low.zero)
               for i in range(n) for j in range(i, n))
    G = identity(low, n)
    for i in range(n):
        for j in range(i):
            value = low.zero
            for k in range(j, i):
                value = low.add(value, low.mul(aa[i][k], G[k][j]))
            G[i][j] = low.scale(value, -1)
    width = 1
    while width < R.e:
        width = min(2*width, R.e)
        rr = R.parent.ring(width)
        G, a = mreduce(rr, G), mreduce(rr, A)
        G = mmul(rr, G, msub(rr, mscale(rr, identity(rr, n), 2), mmul(rr, a, G)))
    return mreduce(R, G)


def contraction_matrix(parent, b, e, operator):
    R = parent.ring(e)
    b = mreduce(R, b)
    if e == 1:
        return b
    m = (e+1)//2
    x = mreduce(R, contraction_matrix(parent, b, m, operator))
    residual = madd(R, msub(R, b, x), operator(x, e))
    divided = [[R.shift_down(v, m) for v in row] for row in residual]
    hi = mreduce(R, contraction_matrix(parent, divided, e-m, operator))
    return madd(R, x, mscale(R, hi, 1 << m))


def balanced_flag(parent, B, gamma, N, P):
    leaves, widths = [], []

    def node(B, offset):
        e = P-offset*N
        R = parent.ring(e)
        B = mreduce(R, B)
        size = len(B)
        widths.append([size, offset, e])
        if size == 1:
            leaves.append((offset, e, B[0][0]))
            return
        k, u = size//2, size-size//2
        b11 = [row[:k] for row in B[:k]]
        b12 = [row[k:] for row in B[:k]]
        b21 = [row[:k] for row in B[k:]]
        b22 = [row[k:] for row in B[k:]]

        def blocks(rr):
            return [mreduce(rr, a) for a in (b11, b12, b21, b22)]

        def Z(rr, y):
            return [[rr.scale(rr.tau(v), pow(gamma, k+i-j, rr.mod))
                     for j, v in enumerate(row)] for i, row in enumerate(y)]

        low = parent.ring(1)
        a11, _, a21, _ = blocks(low)
        Y = mmul(low, a21, minverse(low, a11))
        known = 1
        while known < e:
            target = min(2*known, e)
            rr = parent.ring(target)
            a11, a12, a21, a22 = blocks(rr)
            Y = mreduce(rr, Y)
            z = Z(rr, Y)
            V = madd(rr, a11, mmul(rr, a12, z))
            U = msub(rr, a22, mmul(rr, Y, a12))
            Vinv = minverse(rr, V)
            F = msub(rr, msub(rr, mmul(rr, Y, V), a21), mmul(rr, a22, z))
            rhs = mscale(rr, mmul(rr, F, Vinv), -1)
            C = [[rr.scale(v, pow(gamma, j, rr.mod))
                  for j, v in enumerate(row)] for row in U]
            right = [[rr.scale(v, pow(gamma, k-1-i, rr.mod)) for v in row]
                     for i, row in enumerate(Vinv)]

            def operator(h, width):
                r = parent.ring(width)
                th = [[r.tau(v) for v in row] for row in mreduce(r, h)]
                return mscale(r, mmul(r, mmul(r, mreduce(r, C), th),
                                     mreduce(r, right)), gamma)

            h = contraction_matrix(parent, rhs, target, operator)
            Y = madd(rr, Y, h)
            known = target
        Y = mreduce(R, Y)
        z = Z(R, Y)
        V = madd(R, b11, mmul(R, b12, z))
        U = msub(R, b22, mmul(R, Y, b12))
        F = msub(R, msub(R, mmul(R, Y, V), b21), mmul(R, b22, z))
        assert all(v == R.zero for row in F for v in row)
        node(V, offset)
        node(U, offset+k)

    node(B, 0)
    return leaves, widths


def specialize(acquired, method, parent, L, P):
    if method == "ordinary":
        R = parent.ring(P)
        return [[R.elem(*v) for v in row] for row in acquired['B_coefficients']]
    Q, gamma = acquired['Q'], acquired['gamma']
    RQ, R = parent.ring(Q), parent.ring(P)
    out = []
    for i, row in enumerate(acquired['coefficients']):
        odd = pow(gamma >> 1, i, R.mod)
        out.append([R.scale(R.reduce(RQ.shift_down(RQ.elem(v), i)),
                            pow(odd, -1, R.mod)) for v in row])
    return out


def centered(value, P):
    value %= 1 << P
    return value-(1 << P) if value >= 1 << (P-1) else value


def readout(parent, leaves, gamma, N, L, P, norm_function):
    R = parent.ring(P)
    answer = R.zero
    details = []
    for j, e, alpha in leaves:
        r = parent.ring(e)
        assert r.parent.ring(1).reduce(alpha) == r.parent.ring(1).one
        nr, ni = norm_function([int(c) for c in r.g],
                              [int(c) for c in alpha[0]],
                              [int(c) for c in alpha[1]], e)
        n = r.pow(r.elem([nr], [ni]), N//r.d)
        answer = R.add(answer, R.scale(R.reduce(n), pow(gamma, j*N, R.mod)))
        details.append({'slope': j, 'bits': e, 'norm': [int(nr), int(ni)]})
    assert not answer[1], "Gaussian imaginary trace did not vanish"
    assert answer[0].degree() <= 0
    value = centered((-1)**(L-1)*int(answer[0][0]), P)
    return value, details


def query(acquired, method, N, L, f, a, norm_function, keep_debug=False):
    P = acquired.get('P', acquired.get('output_P', N*(L-1)+2))
    Q = acquired.get('Q', P) if method == 'ah' else P
    start = perf_counter()
    fa, orbit = parameter_polynomial(a, f)
    parent = Presentation(fa, Q)
    setup = perf_counter()
    B = specialize(acquired, method, parent, L, P)
    specialized = perf_counter()
    leaves, widths = balanced_flag(parent, B, acquired['gamma'], N, P)
    flagged = perf_counter()
    value, norms = readout(parent, leaves, acquired['gamma'], N, L, P, norm_function)
    end = perf_counter()
    out = {'a': a, 'value': value, 'degree': len(orbit), 'fa': fa,
           'stages': {'parameter': setup-start, 'specialization': specialized-setup,
                      'flag': flagged-specialized, 'norm': end-flagged,
                      'total': end-start}, 'flag_widths': widths, 'norms': norms}
    if keep_debug:
        out['_debug'] = (parent, B, leaves)
    return out


def ordered_trace(parent, B, gamma, N, L, P):
    """Validation only: independent literal N-step ordered full matrix product."""
    R = parent.ring(P)
    A = [[R.scale(v, pow(gamma, i, R.mod)) for v in row] for i, row in enumerate(B)]
    product = identity(R, L)
    for _ in range(N):
        product = mmul(R, A, product)
        A = [[R.sigma(v) for v in row] for row in A]
    trace = R.zero
    for i in range(L):
        trace = R.add(trace, product[i][i])
    assert not trace[1] and trace[0].degree() <= 0
    return centered((-1)**(L-1)*int(trace[0][0]), P)
