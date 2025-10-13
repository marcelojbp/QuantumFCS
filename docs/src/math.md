# [Mathematical Background](@id math)

As a general scenario, we consider a Lindblad master equation,

```math
\mathcal{L}\rho = \frac{d \rho}{dt} = -i[H, \rho] + \sum_{k=1}^r L_k \, \rho \, L_k^\dagger - \frac{1}{2}\{L^\dagger_k L_k, \rho \}.
```
We introduce $p \le r$ counting fields $N_k$ with weights $\nu_k$. This lets us define the total current,

```math
N(t) = \sum_k \nu_k \, N_k(t)~.
```
Above, $\nu_k$ can be simply integers, e.g. $\pm 1$ to count electrons/photons, or have units to count electric charges ($e$) and  stochastic heat/work (energy). 

We further define the $n$-resolved density matrix $\rho_n(t)$ whose trace equals the probability to have accumulated $n$ jumps at time $t$, $P(n,t) = \operatorname{Tr}[\rho_n(t)]$. Summing over the set of allowed values $\mathcal{N}$ for the total charge $N$, we retrieve the standard density matrix,

```math
\rho(t) = \sum_{n \in \mathcal{N}} \rho_n(t)~.
```

We now consider the Fourier transform of the $n$-resolved density matrix,

```math
\rho_{\chi}(t) = \sum_{n \in \mathcal{N}} e^{i n \chi} \, \rho_n(t)~.
```

$\chi$ is called the counting field and the time evolution of $\rho_{\chi}(t)$ is given by the generalized master equation (GME),

```math
\mathcal{L}_\chi \, \rho_\chi = \bigl(\mathcal{L} + \delta \mathcal{L}_\chi\bigr)\rho_\chi,
```

where

```math
\delta \mathcal{L}_\chi = \sum_{k=1}^p\bigl(1-e^{i\nu_k \chi}\bigr) \, L_k \, (\cdot) \, L_k^\dagger~.
```

**Computing cumulants using recursive methods**

We are ultimately interested in the $n$-th cumulant $\langle\!\langle I^n \rangle\!\rangle$ of the stochastic current,

```math
I(t) = \frac{dN}{dt}~,
```

which we compute through the following recursive scheme,

```math
\langle\!\langle I^n \rangle\!\rangle = \sum_{m=1}^n \binom{n}{m} \, \langle\!\langle \mathbb{1} | \, \mathcal{L}^{(m)} \, | \rho_{\text{ss}}^{(n-m)}(\chi) \rangle\!\rangle~,
```

with the constituents,

```math
| \rho_{\text{ss}}^{(n)}(\chi) \rangle\!\rangle = \mathcal{L}^+ \sum_{m=1}^n \binom{n}{m} \Bigl( \langle\!\langle I^m \rangle\!\rangle - \mathcal{L}^{(m)}\Bigr) | \rho_{\text{ss}}^{(n-m)} \rangle\!\rangle~,
```

```math
\mathcal{L}^{(n)} = \bigl(-i \, \partial_{\chi}\bigr)^n \mathcal{L}_\chi \Big|_{\chi \to 0}~,
```

and $\mathcal{L}^+$ being the Drazin inverse of $\mathcal{L}$.

## References

- **FCS in Lindblad master equations**: [Potts (2024)][potts2024]
- **FCS, recursive scheme, vectorization, Drazin inverse**: [Landi et al. (2024)][landi2024]
- **Detailed exposition of the recursive scheme**: [Flindt et al. (2010)][flindt2010]

[landi2024]: arXiv:2303.04270v4 "Current fluctuations in open quantum systems: Bridging the gap between quantum continuous measurements and full counting statistics"
[potts2024]: https://doi.org/10.48550/arXiv.2406.19206 "Quantum Thermodynamics"
[flindt2010]: arXiv:1002.4506v2 "Counting statistics of transport through Coulomb blockade nanostructures: High-order cumulants and non-Markovian effects"
