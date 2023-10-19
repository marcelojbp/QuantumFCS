"""
    fcscumulants_recursive(H, J, mJ, nC; <keyword arguments>)

Calculate n-th zero-frequency cumulant of full counting statistics using a recursive scheme.

# Arguments
* `H`: Arbitrary operator specifying the Hamiltonian.
* `J`: Vector containing all jump operators which can be of any arbitrary
        operator type.
* `mJ`: Vector containing the monitored jump operators.
* `Jdagger=dagger.(Jdagger)`: Vector containing the hermitian conjugates of the
        jump operators. If they are not given they are calculated automatically.
* `nu`: vector of length 2*length(J) weights for each jump operator. By default down jumps, J, have weight -1 and up-jumps have weight +1.
* `kwargs...`: Further arguments are passed on to the ode solver.
"""
function fcscumulants_recursive(H::AbstractOperator, J, mJ, nC::Int64; rho_ss = steadystate.eigenvector(H, J), nu = vcat(fill(-1, Int(length(J)/2)),fill(1, Int(length(J)/2))))
    l = length(rho_ss)
    IdL = Matrix{ComplexF64}(I, l, l)
    vId = vec(Matrix{ComplexF64}(I, size(rho_ss.data)))'
    Ln = [m_jumps(mJ; n = k, nu = nu) for k=1:nC]
    vrho0 = vec(rho_ss.data)
    vI = 0.0*Vector{Float64}(undef, nC)
    vI[1] = real(vId* Ln[1]*vrho0)
    if nC > 1
        vrho = [Vector{ComplexF64}(undef, l) for j=1:nC]
        vrho[1] = vrho0
        LD = drazin(H, J, rho_ss)
        for n = 2:nC
            vrho[n] = LD*sum(binomial(n-1, m)*(vI[m]*IdL*vrho[n-m] - Ln[m]*vrho[n-m]) for m=1:n-1)
            vI[n] = real(vId*sum(binomial(n, m)*Ln[m]*vrho[n+1-m] for m=1:n))
        end
    end
    return vI
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
    function drazin(H::AbstractOperator, J; ρss = steadystate.eigenvector(H, J))
        d = length(H)
        vId = reshape(Matrix(identityoperator(H).data), d)'
        vss = reshape(Matrix(ρss.data), d)
        vL = Matrix(L.data)
        IdL = Matrix{ComplexF64}(I, d, d)
        Q = IdL - vss*vId
        LD = Q*pinv(vL)*Q
        return LD
    end
 """
    m_jumps(mJ; n=1; nu = vcat(fill(-1, Int(length(J)/2)),fill(1, Int(length(J)/2))))

Calculate the vectorized super-operator ℒ(n) = ∑ₖ (νₖ)ⁿ (Lₖ*)⊗Lₖ.      
# Arguments
* `mJ`: List of monitored jumps
* `n` : Power of the weights νₖ. By default set to 1, since this case appears more often.
* `nu`: Set of weights νₖ, by default set to -1 for emission and +1 for absorption jumps.
"""
    function m_jumps(mJ; n=1, nu = vcat(fill(-1, Int(length(J)/2)),fill(1, Int(length(J)/2))))
        return sum(nu[k]^n*kron(conj(mJ[k].data), mJ[k].data) for k = 1:length(mJ))
    end
    