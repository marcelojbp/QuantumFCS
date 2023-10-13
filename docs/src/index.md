# Documentation 

## Introduction 

As a general scenario, we consider a Lindblad master equations denoted by,

``
\mathcal{L}\rho =\frac{d \rho}{dt} = -i[H, \rho] + \sum_{k=1}^r L_k \rho L_k^\dagger -\frac{1}{2}\left\{L^\dagger_k L_k, \rho \right\}.
``

We introduce ``p \leq r`` counting fields `` N_k `` with weights ``\nu_k ``  which lets us define the total current, 

``
N(t) = \sum_k \nu_k N_k(t) ~ . 
``

We further define the n-resolved density matrix `` \rho_n(t)`` whose trace equals the probability ``P(n,t) = \text{Tr}\left[ \rho_n(t) \right]`` to have accumulated ``n`` jumps at time ``t ``. Summing over the set of allowed values  `` \mathcal{N}`` for the total charge ``N`` , we retrieve the standard density matrix, 

``
\rho = \sum_{n \in \mathcal{N}} \rho_n(t)~.
``

We now consider the Fourier transform of the n-resolved density matrix, 

``
\rho_{\chi}(t) = \sum_{n \in \mathcal{N}} e^{i n \chi} \rho_n(t)~.
``

Here ``\chi`` is the counting field and the dynamics of ``\rho_{\chi}(t)`` is given by the generalized master equation (GME),

``
\mathcal{L}_\chi \rho_\chi = \left(\mathcal{L} + \Delta \mathcal{L}_\chi \right)\rho_\chi,
``

where,

``
\Delta \mathcal{L}_\chi = \sum_{k=1}^p(1-e^{i\nu_k \chi}) L_k \rho L_k^\dagger .
``

## Functions 

```@docs 
fcscumulants_recursive
```


# Examples 

