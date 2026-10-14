@testset "Integer initial data and exact discrete oracles" begin
    p = IVP((u, t)->-u, [1], (0, 1))
    for solver in (RK4(h=0.1), BackwardEuler(h=0.1), RadauIA3(h=0.1))
        sol = solve(p, solver)
        @test sol.t[end] === 1.0
        @test all(isfinite, sol.u[end])
    end
    @test solve(p, RK4(h=0.1)).u[end][1] ≈ exp(-1) rtol=1e-5
    @test solve(p, BackwardEuler(h=0.1)).u[end][1] ≈ (1/1.1)^10 rtol=1e-12
    @test solve(p, RadauIA3(h=0.1)).u[end][1] ≈ exp(-1) rtol=1e-4
    # Backward Euler has first-order global error, not 1% accuracy at h=0.1.
    e1 = abs(solve(p, BackwardEuler(h=0.1)).u[end][1]-exp(-1))
    e2 = abs(solve(p, BackwardEuler(h=0.05)).u[end][1]-exp(-1))
    @test 1.8 < e1/e2 < 2.1
end

@testset "Time precision across all solver families" begin
    for T in (Float32, BigFloat)
        nonlinear = IVP((u, t)->-u, T[1], (zero(T), one(T)))
        linear = IVP(-1.0, T[1], (zero(T), one(T)))
        split = IVP(-1.0, (u, t)->zero.(u), T[1], (zero(T), one(T)))
        for (p, solver) in (
            (nonlinear, RK4(h=0.1)),
            (linear, BackwardEuler(h=0.1)),
            (linear, RadauIA3(h=0.1)),
            (split, IMEXEuler(h=0.1)),
            (linear, ExponentialEuler(h=0.1)),
            (nonlinear, DormandPrince54(h=0.1, εₐ=1e-7, εᵣ=1e-5)),
        )
            sol = solve(p, solver)
            @test eltype(sol.t) === T
            @test eltype(sol.u[end]) === T
            @test sol.t[end] == one(T)
            @test issorted(sol.t)
            @test sol.u[end][1] ≈ exp(-1) rtol=0.05
        end
    end
end
