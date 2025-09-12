```julia
using QuantumFCS

L = #Vectorised Liouvillean
rho_ss = #Steady-state (0-th eigenvalue of the L)

# compute first 3 cumulants
avg, var, skw = fcscumulants_recursive(L, J, 3, rho_ss)