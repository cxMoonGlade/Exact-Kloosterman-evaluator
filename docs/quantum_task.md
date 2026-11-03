# Exact classical simulation of half-period quantum walks on norm graphs

At half their minimum operator period, Laplacian quantum walks on binary norm graphs have nonzero-displacement amplitudes given by hyper-Kloosterman sums. AH evaluates these sums exactly. The resulting probabilities support exact classical sampling of complete vertex labels from any initial vertex.

## Graph, input and output

Let $q=2^N$, $E=\mathbb F_{q^\ell}$ and $n=N\ell$, where $N\ge1$ and $\ell\ge2$. The input supplies monic irreducible polynomials $f$ and $g$, of degrees $N$ and $\ell$, specifying the tower

$$
\mathbb F_q=\mathbb F_2[z]/(f),\qquad E=\mathbb F_q[y]/(g).
$$

A vertex is stored as $\ell$ base-field coefficients, giving an $n$-bit label. Write

$$
Q=q^\ell,\qquad \nu=N_{E/F_q},\qquad d=\frac{Q-1}{q-1}.
$$

The characteristic-two Kollár–Rónyai–Szabó norm graph joins $u,v\in E$ when $\nu(u+v)=1$. It is a cubelike graph of degree $d$. For adjacency matrix $A$, use the Laplacian Hamiltonian $dI-A$ and evolution

$$
U(t)=\exp\!\bigl(-it(dI-A)\bigr).
$$

We study the walk at $t_*=T/2=\pi/2^\ell$, where $T=2\pi/2^\ell$ is its minimum operator period. For an initial vertex $v$, the position distribution is

$$
p_v(u)=\left|\langle u|U(t_*)|v\rangle\right|^2.
$$

A probability query returns $p_v(u)$ as an exact rational number. A sampling call returns the complete $n$-bit label of $u$ with probability $p_v(u)$. Sampling uses independent fair random bits and terminates almost surely.

## Exact amplitudes

Let $tr=Tr_{\mathbb F_q/\mathbb F_2}$. For nonzero $c\in\mathbb F_q$ and rank $r\ge1$, define

$$
Kl_r(c)=
\sum_{\substack{x_1,\ldots,x_r\in\mathbb F_q^\times\\x_1\cdots x_r=c}}
(-1)^{tr(x_1+\cdots+x_r)}.
$$

The arithmetic rank is stored as `L` in the published code and result tables. For the quantum walk, set $r=L=\ell-1$.

A $2$-adic congruence identifies the half-period unitary as

$$
U(t_*)=\mathcal F_E\,
diag_{a\in E}\!\left((-1)^{tr(\nu(a))}\right)
\mathcal F_E,
$$

where

$$
\mathcal F_E|x\rangle=
Q^{-1/2}\sum_{a\in E}(-1)^{Tr_{E/\mathbb F_2}(ax)}|a\rangle.
$$

Writing $b=u+v$, the amplitudes are

$$
\langle u|U(t_*)|v\rangle=
\begin{cases}
(1-d)/Q,&b=0,\\[2pt]
(-1)^r q Kl_r(\nu(b))/Q,&b\ne0.
\end{cases}
$$

Squaring this expression gives each output probability from one exact integer evaluation, with the return probability available directly. For $b\ne0$, the arithmetic parameter is $\nu(b)=Res_y(g,b(y))$, computed in the supplied tower representation.

The Fourier factorization also gives an $O(\ell n^2)$-gate quantum implementation, using $n$ position qubits and polynomially many auxiliary qubits. For $N\ge2$ and $\ell\ge3$, the unitary is non-Clifford. Its position distribution stays far from distributions supported on polynomially many vertices as $n$ grows.

## Exact sampling of complete positions

For $N\ge2$ and $\ell\ge3$, precompute the AH coefficients once for degree $N$ and rank $r$. Deligne's bound $|Kl_r(c)|\le r q^{(r-1)/2}$ gives the squared envelope $B=r^2q^{r-1}$. One sample is generated as follows.

1. Return $v$ with probability $(d-1)^2/Q^2$. Otherwise enter the rejection loop.
2. Propose a uniform nonzero displacement $b\in E^\times$, compute $c=\nu(b)$, and evaluate $Kl_r(c)$ using AH.
3. Accept with probability $Kl_r(c)^2/B$ and return $v+b$. On rejection, repeat step 2.

Each proposal supplies a complete vertex label. Because acceptance depends only on the norm, vertices with the same norm remain equally likely. The return branch and rejection loop together produce exactly $p_v$.

The second-moment identity

$$
\sum_{c\in\mathbb F_q^\times}Kl_r(c)^2
=\frac{1+(q-2)q^r}{q-1}
$$

determines the acceptance rate. If $J$ counts arithmetic evaluations in one complete sampling call, including the initial return branch, then

$$
\mathbb E J=r^2(1-Q^{-1}).
$$

All random decisions use exact integer comparisons. For example, draw a uniform integer below $B$ and compare it with $Kl_r(c)^2$ in step 3.

The boundary cases have direct samplers. At $N=1$, the graph is $K_{2^\ell}$. Return to $v$ has probability $(1-2/Q)^2$, and every other vertex has probability $4/Q^2$. At $\ell=2$, the position distribution is uniform on all $Q$ vertices. Both cases use the explicit probabilities without evaluating AH.

## Uniform total cost

For $N,r\ge2$, AH requires shared preprocessing followed by one evaluation per requested parameter. Each stage has bit cost

$$
\widetilde O\!\left(N^2r^{1+\omega_0}\right),\qquad \omega_0=\log_2 7.
$$

Here $\widetilde O$ suppresses logarithmic factors, and $\omega_0$ is the exponent of Strassen matrix multiplication. One precomputation serves all nonzero parameters at fixed degree and rank. A probability query, including preprocessing and norm computation, therefore takes deterministic bit time

$$
\widetilde O\!\left(N^2\ell^{1+\omega_0}\right).
$$

Generating $S\ge1$ independent complete positions takes expected bit time

$$
\widetilde O\!\left(N^2\ell^{1+\omega_0}\right)
+S\,\widetilde O\!\left(N^2\ell^{3+\omega_0}\right).
$$

The first term covers shared preprocessing; the second includes every rejected proposal, field norm, exact random decision and complete vertex output. The bound is polynomial jointly in $N$ and $\ell$, hence polynomial in the label length $n=N\ell$.

## Listed coherent inputs

The same signed amplitude formula supports explicitly listed coherent inputs with Gaussian-integer coefficients. For $k$ distinct input vertices, a mixture proposal followed by an interference acceptance test uses at most $k[r^2(1-Q^{-1})+k]$ expected AH evaluations in the AH branch. Its expected bit cost is polynomial jointly in $N$, $\ell$, $k$ and the coefficient bit length.

## Published arithmetic implementation

The [core source](../README.md#core-code) implements the integer evaluator; the [benchmarks](../results/README.md) report its arithmetic costs. The complexity bounds above use structured preprocessing and semilinear reduction. The benchmarked implementation computes the Frobenius trace through ordered matrix products.
