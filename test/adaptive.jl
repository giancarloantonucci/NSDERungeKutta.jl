@testset "adaptive stepping (embedded explicit pairs)" begin
    @testset "loud failure without an embedding" begin
        @test_throws ArgumentError ERK(Euler(h=0.1).tableau, 0.1, AP())
        @test_throws ArgumentError DIRK(BackwardEuler(h=0.1).tableau, 0.1, NP(), AP())
        @test_throws ArgumentError IRK(LobattoIIIC2(h=0.1).tableau, 0.1, NP(), AP())
        imex = IMEXEuler(h=0.1)
        @test_throws ArgumentError IERK(imex.implicitableau, imex.explicitableau, 0.1, NP(), AP())
        @test HeunEuler(h=0.1) isa ExplicitRungeKuttaSolver # embedded: constructs fine
    end

    @testset "controller behaviour" begin
        problem = Logistic(0.3, (0.0, 2.0))
        exact = logistic_exact
        loose = solve(problem, DormandPrince54(h=0.1, εᵣ=1e-4))
        tight = solve(problem, DormandPrince54(h=0.1, εᵣ=1e-8))
        err_loose = norm(loose.u[end] - exact(loose.t[end]))
        err_tight = norm(tight.u[end] - exact(tight.t[end]))
        @test err_tight < err_loose
        @test err_tight < 1e-6
        @test numtimesteps(tight) > numtimesteps(loose)

        recording = DormandPrince54(h=0.1, εᵣ=1e-6, save_stepsizes=true)
        rsol = solve(problem, recording)
        hs = rsol.stepsizes
        @test hs !== nothing
        # The record is the REALISED history, not the controller's proposals:
        # accepted[i] is exactly the i-th step taken. (The old code pushed the
        # post-update h; only length checks lived here, so it went unnoticed.)
        @test length(hs.accepted) == numtimesteps(rsol) - 1
        @test hs.accepted ≈ diff(rsol.t)
        # rejected[i]: the h values that FAILED before accepted step i, plus a
        # trailing empty bin opened at the final acceptance.
        @test length(hs.rejected) == length(hs.accepted) + 1
        @test isempty(hs.rejected[end])
        @test all(bin -> all(>(0), bin), hs.rejected)
    end

    @testset "stability-limited vs accuracy-limited stepping" begin
        # Two linear systems engineered to share the SAME exact solution
        # u(t) = [2e⁻ᵗ + sin t, 2e⁻ᵗ + cos t]; Jacobian eigenvalues {-1, -3}
        # versus {-1, -1000}. Equal accuracy therefore demands SIMILAR step
        # counts — every extra step the stiff twin takes is the price of
        # explicit STABILITY, and the rejection count is the stability
        # boundary made visible.
        uexact(t) = [2exp(-t) + sin(t), 2exp(-t) + cos(t)]
        f₁(u, t) = [-2.0 1.0; 1.0 -2.0] * u + [2sin(t), 2(cos(t) - sin(t))]        # {-1, -3}
        f₂(u, t) = [-2.0 1.0; 998.0 -999.0] * u + [2sin(t), 999(cos(t) - sin(t))]  # {-1, -1000}
        function adaptive_run(f)
            solver = Fehlberg45(h=0.1, εᵣ=1e-2, save_stepsizes=true)
            sol = solve(IVP(RHS(f), [2.0, 3.0], (0.0, 10.0)), solver)
            nsteps = numtimesteps(sol) - 1
            nrejections = sum(length, sol.stepsizes.rejected)
            err = maximum(norm(sol.u[n] - uexact(sol.t[n])) for n in eachindex(sol.t))
            return nsteps, nrejections, err
        end
        n₁, r₁, e₁ = adaptive_run(f₁)
        n₂, r₂, e₂ = adaptive_run(f₂)
        @test e₁ < 0.2 && e₂ < 0.2  # both twins accurate: same solution, same tolerance
        @test n₂ > 5n₁              # stability, not accuracy, sets the stiff step count
        @test r₂ > r₁               # the stability boundary shows up as rejections
    end
end

@testset "adaptive storage does not scale with the starting step" begin
    p = Logistic(0.3, (0.0, 100.0))
    solver = DormandPrince54(h=1e-8, εₐ=1e-12, εᵣ=1e-8)
    for dense in (false, true)
        sol = NSDEBase.initialize_solution(p, solver; dense)
        @test length(sol.u) == length(sol.t) == 2
        @test Base.summarysize(sol) < 10_000
        NSDEBase.solve!(sol, p, solver)
        @test sol.t[end] == 100.0
        @test length(sol.t) < 1000
        @test sol.u[end][1] ≈ 1.0 atol=1e-7
        @test !dense || length(sol.k) == length(sol.t)-1
    end
end

@testset "adaptive retry respects local time spacing" begin
    p = Dahlquist([1.0], (6000.0, 6000.01); λ=1e16)
    solver = DormandPrince54(h=eps(6000.0), εₐ=1e-30, εᵣ=1e-30)
    @test_throws StepFailure solve(p, solver)
    # A representable, accepted final step does not need a viable next proposal.
    p = Dahlquist([1.0], (6000.0, nextfloat(6000.0)); λ=0.0)
    sol = solve(p, DormandPrince54(h=eps(6000.0), εₐ=1e-12))
    @test sol.t[end] == p.tspan[2]
    @test sol.u[end] == [1.0]
end

@testset "adaptive initial step may imply an unrepresentable fixed-grid count" begin
    p = Dahlquist([1.0], (0.0, 100.0); λ=0.0)
    solver = DormandPrince54(h=1e-20, εₐ=1e-12)
    @test (p.tspan[2]-p.tspan[1])/solver.stepsize.h > typemax(Int)
    sol = NSDEBase.initialize_solution(p, solver)
    @test length(sol.t) == 2
    NSDEBase.solve!(sol, p, solver)
    @test sol.t[end] == 100.0
    @test sol.u[end] == [1.0]
    @test length(sol.t) < 100
end
