# [Quickstart](@id quickstart)

This package provides tools to compute full-counting-statistics cumulants from a Liouvillian.

## Quickstart examples

```julia
using QuantumFCS

# Build your Liouvillian L and monitored jumps mJ
# L      :: AbstractMatrix{ComplexF64}      (vectorized Liouvillian)
# mJ     :: Vector of sparse jump superoperators
# rho_ss :: Steady-state density matrix (matrix, not vectorized)

# compute first 3 cumulants
c1, c2, c3 = fcscumulants_recursive(L, mJ, 3, rho_ss)

#In case we integrate with QuantumOptics.jl;
using QuantumOptics

# H :: QuantumOptics.Operator Define your Hamiltonian as a Quantum Optics Operator type
# J :: Vector{QuantumOptics.Operator} Define your vector containing the jump opterators 
# mJ:: Vector{QuantumOptics.Operator} Define your vector containing the monitored jump operators 
# rho_ss ::QuantumOptics.Operator steady-state density operator

c1, c2, c3 = fcscumulants_recursive(H, J, mJ, 3, rho_ss)

```

