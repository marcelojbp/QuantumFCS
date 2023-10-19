# Documentation 

## Introduction 

As a general scenario, we consider a Lindblad master equations,

``
\mathcal{L}\rho =\frac{d \rho}{dt} = -i[H, \rho] + \sum_{k=1}^r L_k \rho L_k^\dagger -\frac{1}{2}\left\{L^\dagger_k L_k, \rho \right\}.
``

We introduce ``p \leq r`` counting fields `` N_k `` with weights ``\nu_k ``. This lets us define the total current, 

``
N(t) = \sum_k \nu_k N_k(t) ~ . 
``

We further define the n-resolved density matrix `` \rho_n(t)`` whose trace equals the probability to have accumulated ``n`` jumps at time ``t `` , ``P(n,t) = \text{Tr}\left[ \rho_n(t) \right]``. Summing over the set of allowed values  `` \mathcal{N}`` for the total charge ``N`` , we retrieve the standard density matrix, 

``
\rho = \sum_{n \in \mathcal{N}} \rho_n(t)~.
``

We now consider the Fourier transform of the n-resolved density matrix, 

``
\rho_{\chi}(t) = \sum_{n \in \mathcal{N}} e^{i n \chi} \rho_n(t)~.
``

``\chi`` is called the counting field and the timeevolution of ``\rho_{\chi}(t)`` is given by the generalized master equation (GME),

``
\mathcal{L}_\chi \rho_\chi = \left(\mathcal{L} + \delta \mathcal{L}_\chi \right)\rho_\chi,
``

where,

``
\delta \mathcal{L}_\chi = \sum_{k=1}^p(1-e^{i\nu_k \chi}) L_k \rho L_k^\dagger .
``

**Computing Cumulants using recursive Methods** 

We are ultimately interested in the ``n``-th cumulant ``\langle \langle  I^n \rangle \rangle `` of the stochastic current, 

``
I(t) = \frac{dN}{dt}~,
``


which we compute through the following recursive scheme, 


``
\langle \langle I^n \rangle \rangle = \sum_{m=1}^n \binom{n}{m} \langle \langle \mathbb{1} | \mathcal{L}^{(m)} | \rho_{\text{ss}}^{(n-m)}(\chi) \rangle \rangle ~,
``

with the constituents, 

``
| \rho_{\text{ss}}^{(n)}(\chi) \rangle \rangle = \mathcal{L}^+ \sum_{m=1}^n \binom{n}{m} \left( \langle \langle I^m \rangle \rangle - \mathcal{L}^{(m)}\right) | \rho_{\text{ss}}^{n-m} \rangle \rangle ~,
``

``
\mathcal{L}^{(n)} = \left. \left( -i \partial_{\chi} \right)^n \mathcal{L}_\chi \right|_{\chi \to 0}~,
``

and ``\mathcal{L}^+`` being the Drazin inverse of ``\mathcal{L}``.
## Functions 

```@docs 
fcscumulants_recursive
```

```@docs 
drazin
```

```@docs 
m_jumps
```


# Examples 

