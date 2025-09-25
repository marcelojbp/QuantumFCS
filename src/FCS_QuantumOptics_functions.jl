function fcscumulants_recursive(
    H::Operator,
    J::AbstractVector{<:Operator},
    mJ::AbstractVector{<:Operator},
    nC::Integer,
    rho_ss::AbstractOperator,
    nu::AbstractVector{<:Real},
    )
    L = liouvillian(H, J).data
    return fcscumulants_recursive(L, getfield.(mJ, :data), nC, sparse(rho_ss.data), nu)
end

# Single convenience wrapper (H,J,...). Placed once to avoid method redefinition during precompilation.
function drazin(H::Operator, J, vrho_ss::AbstractVector, vId::AbstractVecOrMat, IdL::AbstractMatrix)
    L = liouvillian(H, J).data
    l = length(vrho_ss)
    IdL_eff = (size(IdL,1) == l && size(IdL,2) == l) ? IdL : Matrix{eltype(L)}(I, l, l)
    vId_row = (size(vId,1) == 1) ? vId : (collect(vec(vId))')
    return drazin(L, vrho_ss, vId_row, IdL_eff)
end

# Compatibility wrapper for tests using QuantumOptics: (H, J, ...) signature
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
    # Delegate to sparse implementation (factorization optional for repeated calls)
    return drazin_apply(L, αs, ρs, vId_vec)
end