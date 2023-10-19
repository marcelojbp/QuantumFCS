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
function fcscumulants_recursive(H::AbstractOperator, J, mJ, nC::Int64; rho_ss = steadystate.eigenvector(H, J), nu = vcat(fill(+1, Int(length(mJ)/2)),fill(-1, Int(length(mJ)/2))))
    l = length(rho_ss)
    # Identity in Liouville space
    IdL = Matrix{ComplexF64}(I, l, l)
    # Vectorired identity operator
    vId = vec(Matrix{ComplexF64}(I, size(rho_ss.data)))'
    # Jump d/dχ n-derivatives, ℒ(n)
    Ln = [m_jumps(mJ; n = k, nu = nu) for k=1:nC]
    # Vectorized steady-state
    vrho0 = vec(rho_ss.data)
    # Empty list which will collect the cumulants
    vI = Vector{Float64}(undef, nC)
    # First cumulant, computed directly
    vI[1] = real(vId* Ln[1]*vrho0)
    # If we are interested in any higher cumulant we start to use the recursive scheme
    if nC > 1
        # Initializing the list of "states" used in the recursion
        vrho = [Vector{ComplexF64}(undef, l) for j=1:nC]
        # We set the first to be the steady-state
        vrho[1] = vrho0
        #  We now can compute the Drazin inverse directly or use the drazin_apply function
        if apply == :false
            # Computing the Drazin inverse
            LD = drazin(H, J, rho_ss)
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
                vrho[n] = drazin_apply(H, J, valpha, rho_ss)
                vI[n] = real(vId*(sum(binomial(n, m)*Ln[m]*vrho[n+1-m] for m=1:n)))
             end 
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
        # Vectorizing the steady-state and the identity
        vId = reshape(Matrix(identityoperator(H).data), d)'
        vss = reshape(Matrix(ρss.data), d)
        # vectorized liouvillian
        vL = matrix(liouvillian(H, J).data)
        # Liouville space's identity
        IdL = Matrix{ComplexF64}(I, d, d)
        # We introduce Q, which projects any vector in the complement of the kernel o L
        Q = IdL - vss*vId
        # The Drazin inverse is computed by projecting the Moore-Penrose pseudo-inverse, computed using pinv.
        LD = Q*pinv(vL)*Q
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

    drazin_apply() 
Calculates the vector resulting from the Drazin inverse being applied to another vector. 

# Arguments 
* `H`: Arbitrary operator specifying the Hamiltonian.
* `J`: Vector containing all jump operators which can be of any arbitrary
    operator type.
* `ρss ` : steady state of the system as a density matrix
* `α` : input state as a density matrix

"""
function drazin_apply(H, J, ρα, ρss = steadystate.master(H,J)[2][1])
    ## Constructing the matrix 

    # constructing the liouvillian from the Hamiltonian and the jump operators 
    L = Matrix(liouvillian(H, J;).data)

    # constructing unitmatrix to append to liouvillian
    unitm = diagm(ones(size(L)[1]))

    # constructing the final matrix consiting of the liouvillian and the unit matrix 
    Mat = cat(L,unitm;dims=1)

    ## Constructing the right hand side 

    # vectorizing the steady state
    ρssvec = vec(Matrix(ρss.data))

    # vectorizing the state that the drazin inverse is being applied to  
    αvec = vec(ρα.data)

    # constructing the right hand side of the linear system 
    rhs = append!(ραvec - ρssvec*sum(ραvec),zeros(size(L)[1]))

    ## returning the result 
    
    return Mat\rhs
end 