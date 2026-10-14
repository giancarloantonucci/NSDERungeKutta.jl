@testset "Nonfinite adaptive trials shrink and retry" begin
    # A stable exact solution, with the RHS undefined far outside its domain.
    # The oversized first trial leaves the domain; smaller trials stay inside.
    calls = Ref(0)
    p = IVP((u, t)->begin
        calls[] += 1
        abs(u[1]) > 10 ? [NaN] : -100 .* u
    end, [1.0], (0.0, 0.1))
    solver = DormandPrince54(h=0.1, εₐ=1e-10, εᵣ=1e-7, save_stepsizes=true)
    sol = @test_logs solve(p, solver)
    @test sol.t[end] == 0.1
    @test sol.u[end][1] ≈ exp(-10) atol=1e-9
    @test first(sol.stepsizes.rejected[1]) == 0.1
    @test sol.stepsizes.rejected[1][2] ≈ 0.02
    # A nonfinite candidate is invalid even if embedded weights cancel exactly.
    testp = IVP((u, t)->u, [1.0], (0.0, 1.0))
    testsolver = HeunEuler(h=0.1)
    work=initialize_cache(testp, testsolver)
    trial=initialize_solution(testp, testsolver)
    NSDERungeKutta.ensure_capacity!(trial, 2, testp.u0)
    trial.u[2] .= Inf
    for stage in work.k
        stage .= 1.0
    end
    next=NSDERungeKutta.adaptivestep!(work, trial, testsolver, testsolver.adaptive, 0.1)
    @test work.n == 1
    @test next ≈ 0.02
    impossible = IVP((u, t)->fill(NaN, length(u)), [1.0], (0.0, 1.0))
    @test_throws StepFailure solve(impossible, DormandPrince54(h=0.1, Mₙ=3))
    # The compensation from each rejected trial must be discarded. Rebuild
    # the time sums using only accepted sizes and compare bitwise.
    t = sol.t[1]
    e=Ref(0.0)
    for (i, h) in enumerate(sol.stepsizes.accepted)
        t = NSDERungeKutta.compensated_sum(t, h, e)
        i == length(sol.stepsizes.accepted) || @test t == sol.t[i + 1]
    end
end

@testset "Complex IVPs infer Jacobians and operator state types" begin
    ps = (IVP(-1im, [1.0], (0, 1)), IVP((u, t)->-1im .* u, [1.0+0im], (0, 1)))
    for p in ps, solver in (RK4(h=0.01), BackwardEuler(h=0.001), RadauIA3(h=0.01))
        sol=solve(p, solver)
        @test eltype(sol.u[end]) <: Complex
        @test sol.u[end][1] ≈ exp(-1im) atol=6e-4
    end
    for stiff in (-1im, (u, t)->-1im .* u)
        p=IVP(stiff, (u, t)->zero.(u), ComplexF64[1], (0, 1))
        sol=solve(p, IMEXEuler(h=0.001))
        @test sol.u[end][1] ≈ exp(-1im) atol=6e-4
    end
    p=IVP((u, t)->-u .^ 3, (u, t)->zero.(u), [20.0], (0.0, 0.001))
    sol=solve(p, IMEXEuler(h=0.001, Mₙ=50, εᵣ=1e-10))
    @test abs(sol.u[end][1]+0.001*sol.u[end][1]^3-20) < 1e-8
    @test_throws ArgumentError NP(εₐ=0, εᵣ=0)
end

@testset "Solution assignment owns copied data" begin
    p=IVP(-1.0, [1.0], (0.0, 1.0))
    solver=RK4(h=0.1)
    sol=solve(p, solver)
    external=[2.0]
    sol[2]=(external, 0.1)
    @test sol.u[2] !== external
    solve!(initialize_cache(p, solver), sol, p, solver)
    @test external == [2.0]
    other=solve(p, solver)
    sol[1]=other[1]
    @test sol.u[1] !== other.u[1]
    sol[1:2]=other[1:2]
    @test sol.u[2] !== other.u[2]
end

@testset "Dense support is explicit" begin
    for F in (RK4, BackwardEuler, LawsonEuler)
        @test_throws ArgumentError solve(Dahlquist(0.5, (0.0, 1.0)), F(h=0.1); dense=true)
    end
    p=IVP((u, t)->[u[2]; -u[1]], [1.0, 0.0], (0.0, 1.0))
    solver=DormandPrince54(h=0.2, εₐ=1e-12, εᵣ=1e-10)
    sol=solve(p, solver; dense=true)
    for t in (0.123, 0.456, 0.789)
        @test sol(t, solver.tableau) ≈ [cos(t), -sin(t)] atol=1e-9
    end
end

@testset "Recipe variables and stability inputs" begin
    sol=solve(IVP((u, t)->-u, [1.0, 2.0, 3.0], (0.0, 1.0)), RK4(h=0.1))
    attrs=Dict{Symbol,Any}(:variables=>[2, 3], :label=>["y" "z"])
    series=RecipesBase.apply_recipe(attrs, sol)
    @test [s.plotattributes[:label] for s in series] == ["y", "z"]
    phase=RecipesBase.apply_recipe(Dict{Symbol,Any}(:variables=>[2, 3]), NSDEBase._PhasePlot(sol))
    @test length(phase[1].args) == 2
    for solver in (IMEXEuler(h=0.1), LawsonEuler(h=0.1))
        for Wrapper in
            (NSDERungeKutta.STABILITY, NSDERungeKutta.STABILITYF, NSDERungeKutta.ORDERSTARF)
            @test_throws ArgumentError RecipesBase.apply_recipe(
                Dict{Symbol,Any}(), Wrapper((solver,))
            )
        end
    end
end
