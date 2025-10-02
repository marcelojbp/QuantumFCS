# [Examples](@id examples)

## Examples

using QuantumOptics
using QuantumFCS

## Simple driven cavity setup

# Define the Hilbert space dimension 
N = 10  # Fock space cutoff

# Create the Fock basis for the cavity mode
basis = FockBasis(N)

# Define operators
a = destroy(basis)  # annihilation operator
a_dag = create(basis)  # creation operator
n = number(basis)   # number operator

# System parameters
ω = 1.0     # cavity frequency
κ = 0.1     # cavity decay rate
Ω = 0.5     # drive strength (coherent drive)

# Hamiltonian (driven cavity)
H = ω * n + Ω * (a_dag + a)  # free evolution + coherent drive

# Jump operators (cavity decay)
J = [sqrt(κ) * a]  # photon loss

# Counting field vector (for one monitored channel)
nu = [-1]

# Monitored jump operators (photons leaving the cavity)
mJ = [sqrt(κ) * a]

# Initial state (vacuum state)
ψ0 = fockstate(basis, 0)
ρ0 = dm(ψ0)

# Calculate steady state
ρss = steadystate.iterative(H, J);

# Calculate the first three cumulants using FCS
c1, c2, c3 = fcscumulants_recursive(H, J, mJ, 3, ρss, nu)

println("\nFull Counting Statistics:")
println("First cumulant (mean photon flux): $c1")
println("Second cumulant (variance): $c2") 
println("Third cumulant (skewness): $c3")