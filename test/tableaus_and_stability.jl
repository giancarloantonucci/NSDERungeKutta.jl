@testset "Butcher tableau sanity" begin
    tableaus = Tuple{String,Any}[]
    for F in [erk_fixed; erk_embedded]
        push!(tableaus, (string(F), F(h=0.1).tableau))
    end
    for F in [dirk_zoo; irk_zoo]
        push!(tableaus, (string(F), F(h=0.1).tableau))
    end
    for F in ierk_zoo
        solver = F(h=0.1)
        push!(tableaus, (string(F) * " (implicit)", solver.implicitableau))
        push!(tableaus, (string(F) * " (explicit)", solver.explicitableau))
    end
    for (name, tab) in tableaus
        @testset "$name" begin
            @test tab.c ≈ vec(sum(tab.A, dims=2)) atol = 1e-12
            @test sum(tab.b) ≈ 1.0 atol = 1e-12
            if !(tab.d isa Nothing)
                @test sum(tab.d) ≈ 1.0 atol = 1e-12
            end
        end
    end
    # The DIRK kernel reads only the lower triangle; a full matrix here would
    # be silently mis-integrated, so pin the shape:
    for F in dirk_zoo
        @test istril(F(h=0.1).tableau.A)
    end
end

@testset "stability function ℛ(z)" begin
    z = -0.5 + 0.3im
    @test stability_function(z, Euler(h=0.1)) ≈ 1 + z
    @test stability_function(z, RungeKutta4(h=0.1)) ≈ sum(z^j / factorial(j) for j in 0:4)
    @test stability_function(z, BackwardEuler(h=0.1)) ≈ 1 / (1 - z)
    @test stability_function(z, ImplicitMidpoint(h=0.1)) ≈ (1 + z / 2) / (1 - z / 2)
    @test stability_function(0.0 + 0.0im, RungeKutta4(h=0.1)) == 1.0

    # A-stability samples: |R| ≤ 1 on the imaginary axis and in the left half-plane
    for solver in (GaussLegendre4(h=0.1), RadauIIA3(h=0.1))
        for y in (0.1, 1.0, 10.0, 100.0)
            @test abs(stability_function(im * y, solver)) ≤ 1 + 1e-10
        end
        for x in (-0.1, -1.0, -10.0), y in (0.0, 1.0, 25.0)
            @test abs(stability_function(x + im * y, solver)) ≤ 1 + 1e-10
        end
    end
    # L-stability behaviour: Radau IIA kills stiff components, Gauss does not
    @test abs(stability_function(-1e6 + 0.0im, RadauIIA3(h=0.1))) < 1e-3
    @test abs(stability_function(-1e6 + 0.0im, GaussLegendre4(h=0.1))) > 0.9

    # Matrix R(Z) against the scalar function on a diagonalisable Z
    P = [1.0 1.0; 0.0 1.0]
    D = [-1.0, -2.0]
    Z = P * Diagonal(D) / P
    for solver in (RungeKutta4(h=0.1), BackwardEuler(h=0.1))
        R_scalar = P * Diagonal([stability_function(complex(λ), solver) for λ in D]) / P
        @test stability_function(Z, solver) ≈ real.(R_scalar) atol = 1e-10
    end
end

@testset "Explicit step-size contract and plot constructors" begin
    for method in
        (Euler, RK4, DormandPrince54, BackwardEuler, GaussLegendre4, IMEXEuler, ExponentialEuler)
        @test_throws UndefKeywordError method()
        solver=method(h=0.1)
        changed=@inferred withstepsize(solver, 0.02)
        @test stepsize(changed)==0.02
        @test stepsize(solver)==0.1
        @test_throws ArgumentError withstepsize(solver, 0)
    end
    @test_throws UndefKeywordError DormandPrince54(εᵣ=1e-8)
    for recipe in (
        NSDERungeKutta.STABILITY,
        NSDERungeKutta.STABILITYF,
        NSDERungeKutta.ORDERSTAR,
        NSDERungeKutta.ORDERSTARF,
    )
        output=RecipesBase.apply_recipe(Dict{Symbol,Any}(), recipe((RK4,)))
        @test !isempty(output)
        @test output[1].args[3](0.0, 0.0)==1.0
    end
end
