
As a general scenario, we consider a Lindblad master equations denoted by,

``
\mathcal{L}\rho =\frac{d \rho}{dt} = -i[H, \rho] + \sum_{k=1}^r L_k \rho L_k^\dagger -\frac{1}{2}\left\{L^\dagger_k L_k, \rho \right\}.
``

From the above equation we introduce $p\leq r$ counting fields with weights $\nu_k$. This defines the generalized master equation (GME),

``
\mathcal{L}_\chi \rho_\chi = \left(\mathcal{L} + \Delta \mathcal{L}_\chi \right)\rho_\chi,
``

where,

``
 \Delta \mathcal{L}_\chi = \sum_{k=1}^p(1-e^{i\nu_k \chi}) L_k \rho L_k^\dagger .
``

Whatever method we use to compute the FCS, we start by vectorizating the GME,

``
\mathcal{L}_\chi\rho_\chi \to   \mathcal{L}_\chi|\rho_\chi \rangle \rangle 
``

```@docs 
fcscumulants_recursive
```