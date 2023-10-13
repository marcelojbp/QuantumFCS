"""
    fcscumulants_recursive(H, J, mJ, nCum; <keyword arguments>)

Calculate n-th zero-frequency cumulant of full counting statistics using a recursive scheme.

# Arguments
* `H`: Arbitrary operator specifying the Hamiltonian.
* `J`: Vector containing all jump operators which can be of any arbitrary
        operator type.
* `mJ`: Vector containing the monitored jumps
* `tol=1e-3`: Tracedistance used as termination criterion.
* `hmin=1e-7`: Minimal time step used in the time evolution.
* `rates=ones(N)`: Vector or matrix specifying the coefficients for the
        jump operators.
* `Jdagger=dagger.(Jdagger)`: Vector containing the hermitian conjugates of the
        jump operators. If they are not given they are calculated automatically.
* `nu`: vector of length 2*length(J) weights for each jump operator. By default down jumps, J, have weight -1 and up-jumps have weight +1.
* `kwargs...`: Further arguments are passed on to the ode solver.
"""
function fcscumulants_recursive(H::AbstractOperator, J, mJ, nCum; tol=1e-3, nu = vcat(-1*ones(length(J), ones(length(J)))))
    return nothing
end
"""
   drazin(H, J; rho_ss = steadystate.eigenvector(H, J))
   
Calculate the Drazin inverse of a Liouvillian defined by the Hamiltonian H and jump operators J.

# Arguments
* `H`: Arbitrary operator specifying the Hamiltonian.
* `J`: Vector containing all jump operators which can be of any arbitrary
        operator type.
* `rho_ss`: Density matrix specifying the steady-state of the Liouvillian. By default, it is found through steadystate.eigenvector. 
        For large matrices the steady-state should be provided, as the best steady-state solver could vary.
"""
function drazin(H::AbstractOperator, J; rho_ss = steadystate.eigenvector(H, J))
        d = length(H)
        vId = reshape(Matrix(identityoperator(H).data), d)'
        vss = reshape(Matrix(rho_ss.data), d)
        vL = Matrix(L.data)
        IdL = Matrix{ComplexF64}(I, d, d)
        Q = IdL - vss*vId
        LD = Q*pinv(vL)*Q
        return LD
    end