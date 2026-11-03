# Exact classical simulation of half-period quantum walks on norm graphs

The task is to compute exact position probabilities and generate complete vertex labels for a continuous-time quantum walk. AH evaluates the variable-rank hyper-Kloosterman sums that determine this distribution.

## Graph, input and output

Let $q=2^N$, $E=\mathbb F_{q^\ell}$ and $n=N\ell$, where $N\ge1$ and $\ell\ge2$. The input gives monic irreducible polynomials specifying the tower

$$
\mathbb F_q=\mathbb F_2[z]/(f),\qquad E=\mathbb F_q[y]/(g),
$$

with degrees $N$ and $\ell$, respectively. A vertex is represented by its $\ell$ base-field coefficients, giving a complete $n$-bit label. Write

$$
Q=q^\ell,\qquad \nu=\operatorname{Norm}_{E/\mathbb F_q},\qquad d=\frac{Q-1}{q-1}.
$$

The characteristic-two Kollár–Rónyai–Szabó norm graph joins $u,v\in E$ when $\nu(u+v)=1$. It is a cubelike graph of degree $d$. For adjacency matrix $A$, use the Laplacian Hamiltonian $dI-A$ and evolution

$$
U(t)=\exp\!\bigl(-it(dI-A)\bigr).
$$

The minimum operator period is $T=2\pi/2^\ell$. The task fixes the evolution time to $t_*=T/2=\pi/2^\ell$. Given an initial vertex $v$, the target law is

$$
p_v(u)=\left|\langle u|U(t_*)|v\rangle\right|^2.
$$

A probability query returns $p_v(u)$ as an exact rational number. A sampling call returns all $n$ bits of a vertex with this probability. Sampling uses independent fair random bits and terminates almost surely.

## From the quantum walk to an arithmetic query

Let $\operatorname{tr}=\operatorname{Tr}_{\mathbb F_q/\mathbb F_2}$. For nonzero $c\in\mathbb F_q$ and rank $r\ge1$, define

$$
\operatorname{Kl}_r(c)=
\sum_{\substack{x_1,\ldots,x_r\in\mathbb F_q^\times\\x_1\cdots x_r=c}}
(-1)^{\operatorname{tr}(x_1+\cdots+x_r)}.
$$

The arithmetic rank is stored as `L` in the published code and result tables. For the quantum walk, set $r=L=\ell-1$.

A $2$-adic congruence identifies the half-period unitary as

$$
U(t_*)=\mathcal F_E\,
\operatorname{diag}_{a\in E}\!\left((-1)^{\operatorname{tr}(\nu(a))}\right)
\mathcal F_E,
$$

where

$$
\mathcal F_E|x\rangle=
Q^{-1/2}\sum_{a\in E}(-1)^{\operatorname{Tr}_{E/\mathbb F_2}(ax)}|a\rangle.
$$

Writing $b=u+v$, the amplitudes are

$$
\langle u|U(t_*)|v\rangle=
\begin{cases}
(1-d)/Q,&b=0,\\[2pt]
(-1)^r q\operatorname{Kl}_r(\nu(b))/Q,&b\ne0.
\end{cases}
$$

Thus one exact integer evaluation determines the probability of any nonzero displacement. Its parameter is computed directly in the supplied tower as $\nu(b)=\operatorname{Res}_y(g,b(y))$.

## Exact sampling of complete positions

For $N\ge2$ and $\ell\ge3$, precompute the AH coefficients once for degree $N$ and rank $r$. Put $B=r^2q^{r-1}$.

1. Return $v$ with probability $(d-1)^2/Q^2$. Otherwise enter the rejection loop.
2. Propose a uniform nonzero displacement $b\in E^\times$, compute $c=\nu(b)$, and evaluate $\operatorname{Kl}_r(c)$ using AH.
3. Accept with probability $\operatorname{Kl}_r(c)^2/B$ and return $v+b$. On rejection, repeat step 2.

Each proposal already supplies a complete vertex label. Acceptance depends only on its norm, so it preserves uniformity within every norm fiber. The resulting law is exactly $p_v$.

Deligne's bound $|\operatorname{Kl}_r(c)|\le r q^{(r-1)/2}$ supplies the envelope. The identity

$$
\sum_{c\in\mathbb F_q^\times}\operatorname{Kl}_r(c)^2
=\frac{1+(q-2)q^r}{q-1}
$$

determines the acceptance rate. If $J$ counts arithmetic evaluations in one complete sampling call, including the initial return branch, then

$$
\mathbb E J=r^2(1-Q^{-1}).
$$

All random decisions use exact integer comparisons. For example, draw a uniform integer below $B$ and compare it with $\operatorname{Kl}_r(c)^2$ in step 3.

The boundary cases use direct distributions. At $N=1$, the graph is $K_{2^\ell}$: return to $v$ has probability $(1-2/Q)^2$, and every other vertex has probability $4/Q^2$. At $\ell=2$, the position distribution is uniform on all $Q$ vertices. These branches require no AH evaluations.

## Uniform total cost

Given the field polynomials, the theoretical AH algorithm has preprocessing and per-evaluation bit costs

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

This includes preprocessing, rejected proposals, field norms, exact random decisions and complete vertex output. The bound is polynomial jointly in $N$ and $\ell$, hence polynomial in the label length $n=N\ell$.

The same signed amplitude formula supports explicitly listed coherent inputs with Gaussian-integer coefficients. For $k$ distinct input vertices, a mixture proposal followed by an interference acceptance test uses at most $k[r^2(1-Q^{-1})+k]$ expected AH evaluations in the AH branch. Its expected bit cost is polynomial jointly in $N$, $\ell$, $k$ and the coefficient bit length.

The quantum evolution also has an $O(\ell n^2)$-gate implementation with $n$ position qubits and polynomially many auxiliary qubits. For $N\ge2$ and $\ell\ge3$, it is non-Clifford and its position law stays far from distributions supported on polynomially many vertices. The exact classical sampler exploits the finite-field arithmetic of this half-period distribution.

## Published arithmetic implementation

The [core source](../README.md#core-code) implements the integer evaluator, and the [saved benchmarks](../results/README.md) measure that arithmetic computation. The theoretical costs above use structured preprocessing and semilinear reduction. The released arithmetic implementation reads out the trace through ordered Frobenius matrix products.
