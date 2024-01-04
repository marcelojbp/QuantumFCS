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
* `nu`: vector of length length(mJ) with weights for each jump operator. By default down jumps, J, have weight +1 and up-jumps have weight -1.
* `apply`: Decide whether we use the drazin_apply or not. By default it is set to :false .
* `kwargs...`: Further arguments are passed on to the ode solver.
"""
function fcscumulants_recursive(H::AbstractOperator, J, mJ, nC::Int64, rho_ss::Operator; apply = :false, nu = vcat(fill(+1, Int(length(mJ)/2)),fill(-1, Int(length(mJ)/2))))
    l = length(rho_ss)
    # Identity in Liouville space
    IdL = Matrix{ComplexF64}(I, l, l)
    # Vectorired identity operator
    vId = vec(Matrix{ComplexF64}(I, size(rho_ss.data)))'
    # Jump d/dχ n-derivatives, ℒ(n)
    Ln = [m_jumps(mJ; n = k, nu = nu) for k=1:nC]
    # Vectorized steady-state
    vrho_ss = vec(rho_ss.data)
    # Empty list which will collect the cumulants
    vI = Vector{Float64}(undef, nC)
    # First cumulant, computed directly
    vI[1] = real(vId* Ln[1]*vrho_ss)
    # If we are interested in any higher cumulant we start to use the recursive scheme
    if nC > 1
        # Initializing the list of "states" used in the recursion
        vrho = [Vector{ComplexF64}(undef, l) for j=1:nC]
        # We set the first to be the steady-state
        vrho[1] = vrho_ss
        #  We now can compute the Drazin inverse directly or use the drazin_apply function
        if apply == :false
            # Computing the Drazin inverse
            LD = drazin(H, J, vrho_ss, vId, IdL)
            for n = 2:nC
            #Computing the "states" 
            valpha = sum(binomial(n-1, m)*(vI[m]*IdL*vrho[n-m] - Ln[m]*vrho[n-m]) for m=1:n-1)
            vrho[n] = LD*valpha
            # and the n-th cumulant
            vI[n] = real(vId*(sum(binomial(n, m)*Ln[m]*vrho[n+1-m] for m=1:n)))
             end
        else
            for n = 2:nC
                # Here we do the same but using the drazin_apply function
                valpha = sum(binomial(n-1, m)*(vI[m]*IdL*vrho[n-m] - Ln[m]*vrho[n-m]) for m=1:n-1)
                vrho[n] = drazin_apply(H, J, valpha, vrho_ss, vId)
                vI[n] = real(vId*(sum(binomial(n, m)*Ln[m]*vrho[n+1-m] for m=1:n)))
             end 
        end
    end
    return vI
    end
"""
    drazin(H, J, vrho_ss, vId, IdL)
    
Calculate the Drazin inverse of a Liouvillian defined by the Hamiltonian H and jump operators J.
 
# Arguments
* `H`: Arbitrary operator specifying the Hamiltonian.
* `J`: Vector containing all jump operators which can be of any arbitrary
         operator type.
* `vrho_ss`: Vectorized density matrix specifying the steady-state of the Liouvillian. 
* `vId` : Vectorized identity operator  
* `IdL` : Identity in Liouville space
 """    
    function drazin(H::AbstractOperator, J, vrho_ss::Vector{ComplexF64}, vId::Adjoint{ComplexF64, Vector{ComplexF64}}, IdL::Matrix{ComplexF64})
        # vectorized liouvillian
        L = Matrix(liouvillian(H, J).data)
        # We introduce Q, which projects any vector in the complement of the kernel o L
        Q = IdL - vrho_ss*vId
        # The Drazin inverse is computed by projecting the Moore-Penrose pseudo-inverse, computed using pinv.
        LD = Q*pinv(L)*Q
        return LD
    end
 """
    m_jumps(mJ; n=1; nu = vcat(fill(-1, Int(length(J)/2)),fill(1, Int(length(J)/2))))

Calculate the vectorized super-operator ℒ(n) = ∑ₖ (νₖ)ⁿ (Lₖ*)⊗Lₖ.      
# Arguments
* `mJ`: List of monitored jumps
* `n` : Power of the weights νₖ. By default set to 1, since this case appears more often.
* `nu`: vector of length length(mJ) with weights for each jump operator. By default down jumps, J, have weight +1 and up-jumps have weight -1.
"""
    function m_jumps(mJ; n=1, nu = vcat(fill(+1, Int(length(mJ)/2)),fill(-1, Int(length(mJ)/2))))
        return sum(nu[k]^n*kron(conj(mJ[k].data), mJ[k].data) for k = 1:length(mJ))
    end
    
""" 

    drazin_apply(H, J, valpha, vrho_ss, vId) 
Calculates the vector resulting from the Drazin inverse being applied to another vector. 

# Arguments 
* `H`: Arbitrary operator specifying the Hamiltonian.
* `J`: Vector containing all jump operators which can be of any arbitrary
    operator type.
* `vrho_ss ` : Steady state of the system as a vectorized density matrix
* `valpha` : Input state as a vectorized density matrix
* `vId`    : Vectorized identity operator

"""
function drazin_apply(H::AbstractOperator, J, valpha::Vector{ComplexF64}, vrho_ss::Vector{ComplexF64}, vId::Adjoint{ComplexF64, Vector{ComplexF64}})
    ## Constructing the matrix 

  # constructing the Liouvillian from the Hamiltonian and the jump operators 
  L = Matrix(liouvillian(H, J;).data)
  
  # constructing the left hand side consiting of the liouvillian and the unit matrix row 
  lhs = cat(L, vId; dims = 1)

  # constructing the right hand side of the linear system 
  rhs = append!(valpha - vrho_ss .* (vId* valpha), 0)

  ## returning the result 
  
  return lhs\rhs

end 