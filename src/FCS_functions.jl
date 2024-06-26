"""
    fcscumulants_recursive(H, J, mJ, nC; <keyword arguments>)

Calculate n-th zero-frequency cumulant of full counting statistics using a recursive scheme.

# Arguments
* `L`: Vectorized Liouvillean matrix
* `mJ`: Vector containing the monitored jump matrices.
* `Jdagger=dagger.(Jdagger)`: Vector containing the hermitian conjugates of the
        jump matrices. If they are not given they are calculated automatically.
* `nu`: vector of length length(mJ) with weights for each jump matrix. By default down jumps, J, have weight +1 and up-jumps have weight -1.
* `apply`: Decide whether we use the drazin_apply or not. By default it is set to :true .
* `iterative`: Further arguments are passed on to the drazin_apply linear solver.
"""
function fcscumulants_recursive(L::SparseMatrixCSC{ComplexF64, Int64}, mJ::Vector, nC::Int64, rho_ss; apply = :true, nu = vcat(fill(+1, Int(length(mJ)/2)),fill(-1, Int(length(mJ)/2))), iterative= :true)
    l = length(rho_ss)
    s = size(rho_ss)
    # Identity in Liouville space
    IdL = spdiagm(ComplexF64.(ones(l)))
    # Vectorired identity operator
    vId = vec(spdiagm(ComplexF64.(ones(s[1]))))'
    # Jump d/dχ n-derivatives, ℒ(n)
    Ln = [m_jumps(mJ; n = k, nu = nu) for k=1:nC]
    # Vectorized steady-state
    vrho_ss = SparseVector(vec(rho_ss))
    # Empty list which will collect the cumulants
    vI = Vector{Float64}(undef, nC)
    # First cumulant, computed directly
    vI[1] = real(vId⋅(Ln[1]*vrho_ss))
    # If we are interested in any higher cumulant we start to use the recursive scheme
    # Initializing the list of "states" used in the recursion
    # We set the first to be the steady-state
    # vrho[1] = vrho_ss
    vrho = [vrho_ss for j=1:nC]
    # vrho[1] = vrho_ss
    #  We now can compute the Drazin inverse directly or use the drazin_apply function
    if apply == :false
        # Computing the Drazin inverse
        LD = drazin(L, vrho_ss, vId, IdL)
        for n = 2:nC
        #Computing the "states" 
        valpha = sum(binomial(n-1, m)*(vI[m]*IdL*vrho[n-m] - Ln[m]*vrho[n-m]) for m=1:n-1)
        vrho[n] = LD*valpha
        # and the n-th cumulant
        append!(vI, real(vId*(sum(binomial(n, m)*Ln[m]*vrho[n+1-m] for m=1:n))))
            end
        else
        for n = 2:nC
            # Here we do the same but using the drazin_apply function
            valpha = sum(binomial(n-1, m)*(vI[m]*IdL*vrho[n-m] - Ln[m]*vrho[n-m]) for m=1:n-1)
            vrho[n] = drazin_apply(L, valpha, vrho_ss, vId; iterative = iterative)
            vI[n] = real(vId⋅(sum(binomial(n, m)*Ln[m]*vrho[n+1-m] for m=1:n)))
            end 
    end

    return vI
end

"""
    fcscumulants_recursive_new(H::Operator, J, mJ, nC, rho_ss::Operator; kargs...)
    Method to integrate with QuantumOptics objets. Instead of the vectorised Liouvillean we pass the Hamiltonian operator, H and a list of jump operators J.
"""
function fcscumulants_recursive(H::Operator, J, mJ, nC, rho_ss::Operator; kargs...)
    L = liouvillian(H, J).data    
    return fcscumulants_recursive(L, getfield.(mJ,:data), nC ,rho_ss.data; kargs...)
end    
"""
    drazin(H, J, vrho_ss, vId, IdL)
    
Calculate the Drazin inverse of a Liouvillian defined by the Hamiltonian H and jump operators J.
 
# Arguments
* `H`: Arbitrary operator specifying the Hamiltonian.
* `J`: Vector containing all jump operators which can be of any arbitrary
         operator type.
* `rho_ss`: Density matrix specifying the steady-state of the Liouvillian. By default, it is found through steadystate.eigenvector. 
         For large matrices the steady-state should be provided, as the best steady-state solver could vary.
 """
function drazin(L, vrho_ss, vId, IdL)
    # vectorized liouvillian
    # We introduce Q, which projects any vector in the complement of the kernel o L
    Q = IdL - vrho_ss*vId
    # The Drazin inverse is computed by projecting the Moore-Penrose pseudo-inverse, computed using pinv.
    LD = Q*pinv(Matrix(L))*Q
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
    return sum(nu[k]^n*kron(conj(mJ[k]), mJ[k]) for k = 1:length(mJ))
end    
"""
    drazin_apply(L::SparseMatrixCSC{ComplexF64, Int64}, alphavec::SparseVector, vrho_ss::SparseVector, vId::AbstractMatrix; iterative=:true)

# Arguments 
* `L`: Vectorized Liouvillean
* `vrho_ss ` : steady state of the system as a vectorized density matrix
* `valpha` : input state as a vectorized density matrix
* `vId`    : vectorized identity operator
* `iterative` : Optional argument to change the linear system solver. As default set to true (Gauss-Seidel), instead of the usual "/" from LinearAlgebra.
"""
function drazin_apply(L::SparseMatrixCSC{ComplexF64, Int64}, alphavec::SparseVector, vrho_ss::SparseVector, vId::AbstractMatrix; iterative=:true)
    # constructing the left hand side consiting of the liouvillian and the unit matrix row 
    lhs = cat(L, vId; dims = 1)
    # constructing the right hand side of the linear system 
    alphavecp = alphavec - vrho_ss .* (vId*alphavec)
    j, v = findnz(alphavecp)
    rhs = SparseVector(length(alphavecp)+1, j, v)
    ## returning the result 
    # return lhs\Vector(rhs)
    # return qr(lhs)\rhs
    if iterative == :true
        # x0, x1 = Krylov.lsmr(lhs, rhs)
        # return  x0
        return IterativeSolvers.gauss_seidel(lhs, Vector(rhs))
    else
         return  lhs\Vector(rhs)
    end
end 
# function drazin_apply(H::AbstractOperator, J, alphavec::Vector{ComplexF64}, vrho_ss::Vector{ComplexF64}, vId::Adjoint{ComplexF64, Vector{ComplexF64}})
#     ## Constructing the matrix 

#   # constructing the Liouvillian from the Hamiltonian and the jump operators 
#   L = Matrix(liouvillian(H, J;).data)
  
#   # constructing the left hand side consiting of the liouvillian and the unit matrix row 
#   lhs = cat(L, vId; dims = 1)

#   # constructing the right hand side of the linear system 
#   rhs = append!(alphavec - vrho_ss .* (vId* alphavec), 0)

#   ## returning the result 
  
#   return lhs\rhs

# end 

# function fcscumulants_recursive(H::AbstractOperator, J, mJ, nC::Int64, rho_ss::Operator; apply = :false, nu = vcat(fill(+1, Int(length(mJ)/2)),fill(-1, Int(length(mJ)/2))))
#     l = length(rho_ss)
#     # Identity in Liouville space
#     IdL = Matrix{ComplexF64}(I, l, l)
#     # Vectorired identity operator
#     vId = vec(Matrix{ComplexF64}(I, size(rho_ss.data)))'
#     # Jump d/dχ n-derivatives, ℒ(n)
#     Ln = [m_jumps(mJ; n = k, nu = nu) for k=1:nC]
#     # Vectorized steady-state
#     vrho_ss = vec(rho_ss.data)
#     # Empty list which will collect the cumulants
#     vI = Vector{Float64}(undef, nC)
#     # First cumulant, computed directly
#     vI[1] = real(vId* Ln[1]*vrho_ss)
#     # If we are interested in any higher cumulant we start to use the recursive scheme
#     if nC > 1
#         # Initializing the list of "states" used in the recursion
#         vrho = [Vector{ComplexF64}(undef, l) for j=1:nC]
#         # We set the first to be the steady-state
#         vrho[1] = vrho_ss
#         #  We now can compute the Drazin inverse directly or use the drazin_apply function
#         if apply == :false
#             # Computing the Drazin inverse
#             LD = drazin(H, J, vrho_ss, vId, IdL)
#             for n = 2:nC
#             #Computing the "states" 
#             valpha = sum(binomial(n-1, m)*(vI[m]*IdL*vrho[n-m] - Ln[m]*vrho[n-m]) for m=1:n-1)
#             vrho[n] = LD*valpha
#             # and the n-th cumulant
#             vI[n] = real(vId*(sum(binomial(n, m)*Ln[m]*vrho[n+1-m] for m=1:n)))
#              end
#         else
#             for n = 2:nC
#                 # Here we do the same but using the drazin_apply function
#                 valpha = sum(binomial(n-1, m)*(vI[m]*IdL*vrho[n-m] - Ln[m]*vrho[n-m]) for m=1:n-1)
#                 vrho[n] = drazin_apply(H, J, valpha, vrho_ss, vId)
#                 vI[n] = real(vId*(sum(binomial(n, m)*Ln[m]*vrho[n+1-m] for m=1:n)))
#              end 
#         end
#     end
#     return vI
#     end

# function drazin(H::AbstractOperator, J, vrho_ss::Vector{ComplexF64}, vId::Adjoint{ComplexF64, Vector{ComplexF64}}, IdL::Matrix{ComplexF64})
#     # vectorized liouvillian
#     L = Matrix(liouvillian(H, J).data)
#     # We introduce Q, which projects any vector in the complement of the kernel o L
#     Q = IdL - vrho_ss*vId
#     # The Drazin inverse is computed by projecting the Moore-Penrose pseudo-inverse, computed using pinv.
#     LD = Q*pinv(L)*Q
#     return LD
# end