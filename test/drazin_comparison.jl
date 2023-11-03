@testset "Drazin comparison" begin
        # Defining the example systems
        basis = FockBasis(4); 
        adagger = create(basis); 
        a = destroy(basis); 
        H = adagger*a; 
        J = [a] ; 

        # Defining the input state and calculating the steady state 
        α = Complex((1+im)/sqrt(2)); 
        ρα = coherentstate(basis, α) ⊗ dagger(coherentstate(basis,α));
        valpha = vec(Matrix(ρα.data));
        rho_ss = steadystate.master(H, J)[2][2];
        vrho_ss = vec(Matrix(rho_ss.data));

        # Defining the vectorized identity (this will usually be passed from the higher level function)
        vId = vec(Matrix{ComplexF64}(I, size(rho_ss.data)))'

        # testing whether drazin_apply and drazin yield the same result
        @test norm(drazin_apply(H, J, valpha, vrho_ss, vId)-drazin(H, J; rho_ss)*valpha) ≈ 0 atol=10^(-10)
end
