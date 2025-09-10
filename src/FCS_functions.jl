"""
    fcscumulants_recursive(L, mJ, nC; <keyword arguments>)

Calculate n-th zero-frequency cumulant of full counting statistics using a recursive scheme.

# Arguments
* `L`: Vectorized Liouvillian matrix (sparse, ComplexF64)
* `mJ`: Vector containing the monitored jump matrices (sparse operators in vectorized representation).
* `nC`: Number of cumulants to be calculated.
* `nu`: Vector of length `length(mJ)` with weights for each jump matrix. By default down jumps have weight -1 and up-jumps have weight +1.
"""
function fcscumulants_recursive(
    L::SparseMatrixCSC{ComplexF64, Int},
    mJ::AbstractVector{<:SparseMatrixCSC{ComplexF64, Int}},
    nC::Integer,
    rho_ss::SparseMatrixCSC{ComplexF64, Int};
    nu = vcat(fill(-1, Int(length(mJ) ÷ 2)), fill(+1, Int(length(mJ) ÷ 2))),
)
    # Dimensions
    n = size(rho_ss, 1)           # matrix side
    l = n * n                     # vectorized length

    # Cached factorization of Liouvillian
    F = lu(L)

    # Vectorized identity (diagonal entries of an n×n identity under vec)
    # Indices: 1:(n+1):l in column-major vectorization
    diag_idx = collect(1:(n+1):l)
    vId = SparseVector{ComplexF64,Int}(l, diag_idx, fill(1.0 + 0.0im, n))

    # d/dχ n-derivatives ℒ(n)
    Ln = [m_jumps(mJ; n = k, nu = nu) for k = 1:nC]

    # Vectorized steady state, normalized
    vrho_ss = SparseVector(vec(rho_ss ./ tr(rho_ss)))

    # Outputs
    vI = Vector{Float64}(undef, nC)

    # First cumulant: I₁ = Re( vId⋅(ℒ(1)*ρ_ss) )
    vI[1] = real(dot(vId, Ln[1] * vrho_ss))

    # States used in recursion
    vrho = Vector{SparseVector{ComplexF64, Int}}(undef, nC)
    vrho[1] = vrho_ss

    # --- Work buffers (reused) ---
    # dense buffers of length l to avoid repeated allocations
    tmp = Vector{ComplexF64}(undef, l)        # for Ln[m] * vrho[·]
    αbuf = zeros(ComplexF64, l)               # accumulates valpha densely

    # main recursion
    for ncur = 2:nC
        # Build valpha = Σ_{m=1}^{n-1} binom(n-1,m) * ( vI[m]*vrho[n-m] - Ln[m]*vrho[n-m] )
        fill!(αbuf, 0)
        for m = 1:(ncur-1)
            c = binomial(ncur - 1, m)
            # αbuf += c * vI[m] * vrho[n-m]
            sv = vrho[ncur - m]
            @inbounds for k in 1:nnz(sv)
                i = rowvals(sv)[k]
                αbuf[i] += c * vI[m] * nonzeros(sv)[k]
            end
            # tmp = Ln[m] * vrho[n-m]; αbuf -= c * tmp
            mul!(tmp, Ln[m], sv)             # SparseMatrix * SparseVector -> dense tmp
            @inbounds @simd for i in eachindex(tmp)
                αbuf[i] -= c * tmp[i]
            end
        end

        # y_n = L^D * valpha  (project+gauge inside drazin_apply), sparse result
        vrho[ncur] = drazin_apply(L, αbuf, vrho_ss, vId; F = F, rtol = 1e-12)

        # I_n = Re( Σ_{m=1}^n binom(n,m) * vId⋅(Ln[m] * vrho[n+1-m]) )
        acc = 0.0
        for m = 1:ncur
            mul!(tmp, Ln[m], vrho[ncur + 1 - m])
            acc += binomial(ncur, m) * real(dot(vId, tmp))
        end
        vI[ncur] = acc
    end

    return vI
end

function fcscumulants_recursive(
    L::Matrix{ComplexF64},
    mJ::AbstractVector{<:SparseMatrixCSC{ComplexF64,Int}},
    nC::Integer,
    rho_ss::Union{SparseMatrixCSC{ComplexF64,Int}, Matrix{ComplexF64}};
    nu = vcat(fill(+1, Int(length(mJ) ÷ 2)), fill(-1, Int(length(mJ) ÷ 2))),
)
    # Dimensions
    n = size(rho_ss, 1)
    l = n*n

    # Factorize dense L once (pivoted LU)
    F = lu(L)

    # Vectorized identity as a sparse vector: indices 1:(n+1):l (column-major)
    diag_idx = collect(1:(n+1):l)
    vId = SparseVector{ComplexF64,Int}(l, diag_idx, fill(1.0 + 0.0im, n))

    # ℒ(n) derivative matrices (still sparse is fine)
    Ln = [m_jumps(mJ; n = k, nu = nu) for k = 1:nC]

    # Vectorized steady-state, normalized, **dense** state vector
    trρ = tr(rho_ss)
    vrho1_dense = vec(Matrix(rho_ss) ./ trρ)  # Vector{ComplexF64}

    # Output cumulants
    vI = Vector{Float64}(undef, nC)

    # I₁ = Re( vId ⋅ (ℒ(1) * ρ_ss) )
    tmp = Vector{ComplexF64}(undef, l)             # dense work buffer length l
    mul!(tmp, Ln[1], vrho1_dense)                  # SparseMatrix * dense -> dense
    vI[1] = real(dot(vId, tmp))

    # States used in recursion, all **dense** for type stability
    vρ = Vector{Vector{ComplexF64}}(undef, nC)
    vρ[1] = vrho1_dense

    # Another dense buffer for α accumulation
    αbuf = zeros(ComplexF64, l)

    for ncur = 2:nC
        # valpha = Σ_{m=1}^{n-1} C(n-1,m) * ( vI[m]*vρ[n-m] - Ln[m]*vρ[n-m] )
        fill!(αbuf, 0)
        for m = 1:(ncur-1)
            c = binomial(ncur - 1, m)

            # αbuf += c * vI[m] * vρ[n-m]
            vnm = vρ[ncur - m]
            @inbounds @simd for i in eachindex(vnm)
                αbuf[i] += c * vI[m] * vnm[i]
            end

            # tmp = Ln[m] * vρ[n-m]; αbuf -= c * tmp
            mul!(tmp, Ln[m], vnm)
            @inbounds @simd for i in eachindex(tmp)
                αbuf[i] -= c * tmp[i]
            end
        end

        # y_n = L^D * valpha  (project+gauge inside drazin_apply)
        vρ[ncur] = drazin_apply(L, αbuf, SparseVector(vrho1_dense), vId; F = F)

        # I_n = Re( Σ_{m=1}^n C(n,m) * vId⋅(Ln[m] * vρ[n+1-m]) )
        acc = 0.0
        for m = 1:ncur
            mul!(tmp, Ln[m], vρ[ncur + 1 - m])
            acc += binomial(ncur, m) * real(dot(vId, tmp))
        end
        vI[ncur] = acc
    end

    return vI
end

function fcscumulants_recursive(
    H::Operator,
    J::AbstractVector{<:Operator},
    mJ::AbstractVector{<:Operator},
    nC::Integer,
    rho_ss::AbstractOperator; kargs...)
    L = liouvillian(H, J).data    
    return fcscumulants_recursive(L, getfield.(mJ, :data), nC, rho_ss.data; kargs...)
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
 
# Convenience wrapper to keep compatibility with tests expecting (H, J, ...) signature
function drazin(H::Operator, J, vrho_ss::AbstractVector, vId::AbstractVecOrMat, IdL::AbstractMatrix)
    L = liouvillian(H, J).data
    vId_vec = collect(vec(vId))
    return drazin(L, vrho_ss, vId_vec, IdL)
end
 """
    m_jumps(mJ; n=1; nu = vcat(fill(-1, Int(length(J)/2)),fill(1, Int(length(J)/2))))

Calculate the vectorized super-operator ℒ(n) = ∑ₖ (νₖ)ⁿ (Lₖ*)⊗Lₖ.      
# Arguments
* `mJ`: List of monitored jumps
* `n` : Power of the weights νₖ. By default set to 1, since this case appears more often.
* `nu`: vector of length length(mJ) with weights for each jump operator. By default down jumps, J, have weight +1 and up-jumps have weight -1.
"""
function m_jumps(mJ::AbstractVector{<:SparseMatrixCSC{ComplexF64, Int}}; n::Integer = 1, nu = vcat(fill(+1, Int(length(mJ) ÷ 2)), fill(-1, Int(length(mJ) ÷ 2))))
    # Sum of sparse Kronecker products stays sparse; element types remain ComplexF64
    return sum(nu[k]^n * kron(conj(mJ[k]), mJ[k]) for k = 1:length(mJ))
end    
# """
#     drazin_apply(L::SparseMatrixCSC{ComplexF64, Int}, alphavec::SparseVector{ComplexF64, Int}, vrho_ss::SparseVector{ComplexF64, Int}, vId::AbstractVector{ComplexF64}; iterative=false, tol=1e-8, maxiter=10_000)

# # Arguments 
# * `L`: Vectorized Liouvillean
# * `vrho_ss ` : steady state of the system as a vectorized density matrix
# * `valpha` : input state as a vectorized density matrix
# * `vId`    : vectorized identity operator (as a vector)
# * `iterative` : Optional argument to change the linear system solver. If true, use GMRES from IterativeSolvers; otherwise use the `\` operator.
# """
# function drazin_apply(
#     L::SparseMatrixCSC{ComplexF64, Int},
#     alphavec::SparseVector{ComplexF64, Int},
#     vrho_ss::SparseVector{ComplexF64, Int},
#     vId::AbstractVector{ComplexF64};
#     # iterative::Bool = false,
#     tol::Real = 1e-8,
#     maxiter::Integer = 10_000,
# )
#     # Project RHS onto the range of L: α' = α - ρ_ss (vId ⋅ α)
#     αp = alphavec - vrho_ss .* (dot(vId, alphavec))
#     # y = L \ αp
#     y = sparse(L \ Vector(αp))
#     # Enforce the gauge (normalization) constraint vId ⋅ y = 0
#     return y - vrho_ss .* (dot(vId, y))
# end 

"""
    drazin_apply(L, α, ρ, vId; F=nothing, rtol=1e-12, atol=0.0)

Apply the (projected) Drazin inverse of `L` to `α` by solving a linear system.
- `L::SparseMatrixCSC{ComplexF64,Int}` (can reuse with `F = factorize(L)`)
- `α::SparseVector{ComplexF64,Int}`
- `ρ::SparseVector{ComplexF64,Int}`   (the steady-state vector)
- `vId::AbstractVector{ComplexF64}`   (gauge vector)
- `F::Union{Nothing,Factorization}` (optional factorization of `L` to reuse)
Returns a `SparseVector{ComplexF64,Int}` (type-stable).
"""
function drazin_apply(L::SparseMatrixCSC{ComplexF64,Int},
                      α::SparseVector{ComplexF64,Int},
                      ρ::SparseVector{ComplexF64,Int},
                      vId::AbstractVector{ComplexF64};
                      F::Union{Nothing,Factorization}=nothing,
                      rtol::Float64=1e-12,
                      atol::Float64=0.0)

    # Project RHS onto range(L): α' = α - ρ * (vId⋅α)
    sα = dot(vId, α)                      # Complex scalar (uses conjugate in dot)
    αp = α - (ρ .* sα)                    # stays SparseVector

    # y = L \ αp   (reuse factorization if provided)
    y = F === nothing ? (L \ Vector(αp)) : (F \ Vector(αp))  # dense Vector{ComplexF64}

    # Enforce gauge: y ← y - ρ * (vId⋅y), but do it *in-place* without densifying ρ
    sy = dot(vId, y)
    @inbounds for k in 1:nnz(ρ)
        i = rowvals(ρ)[k]                 # index of k-th stored entry
        y[i] -= sy * nonzeros(ρ)[k]       # y[i] -= sy * ρ[i]
    end

    # One-pass sparsification with absolute/relative threshold
    thr = max(atol, rtol * norm(y, Inf))
    nzI = Int[]; nzV = ComplexF64[]
    @inbounds for i in eachindex(y)
        yi = y[i]
        if abs(yi) > thr
            push!(nzI, i); push!(nzV, yi)
        end
    end
    return sparsevec(nzI, nzV, length(y))  # SparseVector{ComplexF64,Int}
end


# Add a dense-RHS method
function drazin_apply(L::SparseMatrixCSC{ComplexF64,Int},
                      α::AbstractVector{ComplexF64},       # DENSE RHS here
                      ρ::SparseVector{ComplexF64,Int},
                      vId::SparseVector{ComplexF64,Int};
                      F::Union{Nothing,Factorization}=nothing,
                      rtol::Float64=1e-12,
                      atol::Float64=0.0)

    # Project RHS onto range(L): α' = α - ρ * (vId⋅α)
    sα = dot(vId, α)
    # work: copy α into y (we’ll reuse it as the solve output as well)
    y = copy(α)
    @inbounds for k in 1:nnz(ρ)
        i = rowvals(ρ)[k]
        y[i] -= sα * nonzeros(ρ)[k]
    end

    # Solve once (reuse factorization if given)
    if F === nothing
        y = L \ y
    else
        y = F \ y
    end

    # Enforce gauge: y ← y - ρ * (vId⋅y)
    sy = dot(vId, y)
    @inbounds for k in 1:nnz(ρ)
        i = rowvals(ρ)[k]
        y[i] -= sy * nonzeros(ρ)[k]
    end

    # One-pass sparsification
    thr = max(atol, rtol * norm(y, Inf))
    nzI = Int[]; nzV = ComplexF64[]
    @inbounds for i in eachindex(y)
        yi = y[i]
        if abs(yi) > thr
            push!(nzI, i); push!(nzV, yi)
        end
    end
    return sparsevec(nzI, nzV, length(y))
end


# vId provided as a 1×N row (e.g., SparseMatrixCSC)
function drazin_apply(L::SparseMatrixCSC{T,Int},
                      x::AbstractVector{T},
                      vrho_ss::AbstractVector{T},
                      vId_row::AbstractMatrix{T}) where {T<:Number}
    @assert size(vId_row,1) == 1
    vId_vec = vec(permutedims(vId_row)) :: Vector{T}
    return drazin_apply(L, x, vrho_ss, vId_vec)
end


"""
    drazin_apply(L, α, ρ, vId; F=nothing)

Fast apply of the (projected) Drazin inverse for **dense** L.
- L :: Matrix{ComplexF64}
- α :: AbstractVector{ComplexF64}         (dense RHS)
- ρ :: SparseVector{ComplexF64,Int}       (steady state, sparse)
- vId :: SparseVector{ComplexF64,Int}     (vectorized identity, sparse)
- F :: Union{Nothing,LU}                  (cached `lu(L)`)

Returns Vector{ComplexF64} (type-stable, dense).
"""
function drazin_apply(L::Matrix{ComplexF64},
                      α::AbstractVector{ComplexF64},
                      ρ::SparseVector{ComplexF64,Int},
                      vId::SparseVector{ComplexF64,Int};
                      F::Union{Nothing,LU}=nothing)

    @assert size(L,1) == size(L,2) == length(α) == length(vId)

    # α' = α - ρ * (vId⋅α)   (project RHS onto range(L))
    sα = dot(vId, α)
    y = copy(α)
    @inbounds for k = 1:nnz(ρ)
        i = rowvals(ρ)[k]
        y[i] -= sα * nonzeros(ρ)[k]
    end

    # Solve Ly = α'  (reuse factorization if provided)
    y = F === nothing ? (L \ y) : (F \ y)

    # Gauge: y ← y - ρ * (vId⋅y)
    sy = dot(vId, y)
    @inbounds for k = 1:nnz(ρ)
        i = rowvals(ρ)[k]
        y[i] -= sy * nonzeros(ρ)[k]
    end

    return y  # Vector{ComplexF64}
end

# Compatibility wrapper for tests: (H, J, ...) signature
function drazin_apply(
    H::Operator,
    J,
    alphavec::AbstractVector,
    vrho_ss::AbstractVector,
    vId::AbstractVecOrMat;
    # iterative::Bool = false,
    tol::Real = 1e-8,
    maxiter::Integer = 10_000,
)
    L = liouvillian(H, J).data
    αs = SparseVector(alphavec)
    ρs = SparseVector(vrho_ss)
    vId_vec = collect(vec(vId))
    return drazin_apply(L, αs, ρs, vId_vec;  tol = tol, maxiter = maxiter)
end


# Old and simpler implementations

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



# function drazin(H::AbstractOperator, J, vrho_ss::Vector{ComplexF64}, vId::Adjoint{ComplexF64, Vector{ComplexF64}}, IdL::Matrix{ComplexF64})
#     # vectorized liouvillian
#     L = Matrix(liouvillian(H, J).data)
#     # We introduce Q, which projects any vector in the complement of the kernel o L
#     Q = IdL - vrho_ss*vId
#     # The Drazin inverse is computed by projecting the Moore-Penrose pseudo-inverse, computed using pinv.
#     LD = Q*pinv(L)*Q
#     return LD
# end


# function fcscumulants_recursive(
#     L::SparseMatrixCSC{ComplexF64, Int},
#     mJ::AbstractVector{<:SparseMatrixCSC{ComplexF64, Int}},
#     nC::Integer,
#     rho_ss::SparseMatrixCSC{ComplexF64, Int};
#     # rho_ss::Matrix{ComplexF64};
#     nu = vcat(fill(+1, Int(length(mJ) ÷ 2)), fill(-1, Int(length(mJ) ÷ 2))),
#     # iterative::Bool = false,
#     )
#     l = length(rho_ss)
#     s = size(rho_ss)
#     # Caches the LU-factorisation of the Liouvillean
#     F = lu(L)
#     # Identity in Liouville space
#     IdL = spdiagm(ComplexF64.(ones(l)))
#     # Vectorired identity operator
#     # Keep as a vector for type-stable dot products and to avoid dense row materialization
#     vId = SparseVector(vec(spdiagm(ComplexF64.(ones(s[1])))))
#     # Jump d/dχ n-derivatives, ℒ(n)
#     Ln = [m_jumps(mJ; n = k, nu = nu) for k=1:nC]
#     # Vectorized steady-state
#     # As in Flindt's paper, we enforce normalisation of the steady state
#     vrho_ss = SparseVector(vec(rho_ss./tr(rho_ss)))
#     # Empty list which will collect the cumulants
#     vI = Vector{Float64}(undef, nC)
#     # First cumulant, computed directly
#     vI[1] = real(dot(vId, Ln[1] * vrho_ss))
#     # If we are interested in any higher cumulant we start to use the recursive scheme
#     # Initializing the list of "states" used in the recursion
#     # We set the first to be the steady-state
#     vrho = Vector{SparseVector{ComplexF64, Int}}(undef, nC)
#     vrho[1] = vrho_ss
#             for n = 2:nC
#                 # Here we do the same but using the drazin_apply function
#                 valpha = sum(binomial(n-1, m) * (vI[m] * vrho[n-m] - Ln[m] * vrho[n-m]) for m = 1:n-1)
#                 vrho[n] = drazin_apply(L, valpha, vrho_ss, vId; F = F)
#                 vI[n] = real(sum(binomial(n, m) * dot(vId, Ln[m] * vrho[n+1-m]) for m = 1:n))
#             end 
#         # end

#     return vI
# end