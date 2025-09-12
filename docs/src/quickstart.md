# [Quickstart](@id quickstart)

This package provides tools to compute full-counting-statistics cumulants from a Liouvillian.

Example skeleton (replace with your model):

```julia
using QuantumFCS

# Build your Liouvillian L and monitored jumps mJ
# L      :: AbstractMatrix{ComplexF64}      (vectorized Liouvillian)
# mJ     :: Vector of sparse jump superoperators
# rho_ss :: Steady-state density matrix (matrix, not vectorized)

# compute first 3 cumulants
c1, c2, c3 = fcscumulants_recursive(L, mJ, 3, rho_ss)
```
