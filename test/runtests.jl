using QuantumFCS
using QuantumOptics
using Test
using LinearAlgebra

@testset "QuantumFCS.jl" begin
    include("drazin_comparison.jl")
    include("drazin_inverse.jl")
end

