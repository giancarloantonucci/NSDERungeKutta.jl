# Heun-Euler on u'=u has embedded difference exactly h²/2 per unit start.
function heun_trial(h; dim=1, atol=0.01, rtol=0.0)
    p = IVP((u, t)->u, ones(dim), (0.0, max(1.0, h)))
    solver = HeunEuler(; h=h, εₐ=atol, εᵣ=rtol)
    cache = initialize_cache(p, solver)
    sol = initialize_solution(p, solver)
    NSDERungeKutta.ensure_capacity!(sol, 2, p.u0)
    NSDERungeKutta.step!(cache, sol, p, solver, h)
    proposal = NSDERungeKutta.adaptivestep!(cache, sol, solver, solver.adaptive, h)
    return cache.n == 2, proposal
end
@testset "Embedded local error and component scaling" begin
    accepted, proposal = heun_trial(0.05)
    @test accepted # h²/2 = 0.00125 < 0.01
    @test proposal ≈ 0.1 # maximum growth
    accepted, proposal = heun_trial(2.0; atol=1.5)
    @test !accepted # error 2 exceeds tolerance 1.5 even for h > 1
    @test proposal ≈ 2*0.9/sqrt(2/1.5)
    @test first(heun_trial(1.0; atol=0.5)) # equality is accepted
    for d in (1, 1000)
        @test !first(heun_trial(0.7; dim=d, atol=0.0, rtol=0.1))
        @test first(heun_trial(0.7; dim=d, atol=0.0, rtol=0.14)) # new state sets scale
        @test last(heun_trial(0.7; dim=d, atol=0.0, rtol=0.1)) ≈ 0.7*0.9/sqrt(0.245/0.1945) rtol=1e-12
    end
end
@testset "Adaptive zero starts, accuracy and dimension independence" begin
    zero_problem = IVP((u, t)->zero.(u), [0.0], (0.0, 1.0))
    zero_sol = @test_logs solve(zero_problem, DormandPrince54(h=0.1, εₐ=0.0, εᵣ=1e-8))
    @test zero_sol.u[end] == [0.0]
    @test length(zero_sol.t) <= 6
    p = IVP(LRHS(hcat(-1.0), t->[sin(t)]), [0.0], (0.0, 1.0))
    sol = @test_logs solve(p, DormandPrince54(h=0.1, εₐ=1e-10, εᵣ=1e-8))
    @test sol.u[end][1] ≈ (sin(1)-cos(1)+exp(-1))/2 atol=1e-9
    histories = []
    for dim in (1, 1000)
        q = IVP((u, t)->-u, ones(dim), (0.0, 1.0))
        result = solve(q, DormandPrince54(h=0.1, εₐ=1e-12, εᵣ=1e-8))
        @test maximum(abs.(result.u[end] .- exp(-1))) < 1e-8
        push!(histories, result.t)
    end
    @test histories[1] ≈ histories[2] rtol=1e-6
    for (a, b, m) in
        ((0.0, 0.0, 1), (-1.0, 1.0, 1), (1.0, -1.0, 1), (NaN, 1.0, 1), (1.0, Inf, 1), (1.0, 1.0, 0))
        @test_throws ArgumentError NSDERungeKutta.AdaptiveParameters(εₐ=a, εᵣ=b, Mₙ=m)
    end
    @test DormandPrince54(h=0.1).adaptive.εₐ == 1e-12
end
@testset "Getting started documentation example" begin
    problem = IVP((u, t)->[u[2]; -sin(u[1])], [0.0, π/4], (0.0, 10.0))
    for solver in (RK4(h=1e-2), DormandPrince54(h=1e-2, εᵣ=1e-8))
        solution = solve(problem, solver)
        @test solution.t[end] == 10.0
        @test all(isfinite, solution(5.0))
    end
end
