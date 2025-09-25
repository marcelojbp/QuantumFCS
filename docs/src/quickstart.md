# [Quickstart](@id quickstart)

This package provides tools to compute full-counting-statistics cumulants from a Liouvillian.

## Installation
To install the package, in the Julia REPL, 
```julia
] add https://github.com/marcelojbp/QuantumFCS
```

## Quickstart examples

```julia
using QuantumFCS

# Build your Liouvillian L and monitored jumps mJ
# L     : Complex sparse/dense matrix   (vectorized Liouvillian)
# mJ    : Vector of sparse jump super-operators you want to monitor
# nC    : Number of cumulats to be computed
# rho_ss : Steady-state density matrix (matrix, not vectorized)
# nu     : Vector of weighs (same length as mJ) for the monitored jumps

mJ = [sqrt(kappa) * a, sqrt(kappa) * a_dagger] #Monitoring loss and injection of photons
nu = [-1, 1] #We count -1 if we anihilate and +1 if we create
# compute first 3 cumulants
c1, c2, c3 = fcscumulants_recursive(L, mJ, 3, rho_ss, nu)

#In case you want to use QuantumOptics.jl;
using QuantumOptics

# H :: QuantumOptics.Operator Define your Hamiltonian as a Quantum Optics Operator type
# J :: Vector{QuantumOptics.Operator} Define your vector containing the jump opterators 
# mJ:: Vector{QuantumOptics.Operator} Define your vector containing the monitored jump operators 
# rho_ss ::QuantumOptics.Operator steady-state density operator

mJ = [sqrt(kappa) * a, sqrt(kappa) * a_dagger] #Same as above, but here a and a_dagger are QuantumOptics.Operators
nu = [-1, 1] #Same as above

c1, c2, c3 = fcscumulants_recursive(H, J, mJ, 3, rho_ss, nu)

```

