@testset "a dividing step ends exactly on tN" begin
    # 0.1 is not representable, so three steps of 0.1 from 0 need not land on
    # 0.3 bit for bit; a landing within rounding is snapped onto tN so the
    # solve neither takes a spurious fourth step nor stops an ulp past the end.
    for (tspan, h, nsteps) in (
        ((0.0, 0.3), 0.1, 3),
        ((0.0, 1.0), 1e-3, 1000),
        ((0.25, 0.5), 1e-3, 250),
        ((0.0, 0.7), 0.1, 7),
    )
        sol = solve(Logistic(0.3, tspan), RK4(h=h))
        @test sol.t[end] == tspan[2]
        @test length(sol.t) == nsteps + 1
    end
    over = solve(Logistic(0.3, (0.0, 0.3)), RK4(h=0.2))     # 0.2 does not divide 0.3
    @test over.t[end] == 0.3 && length(over.t) == 3          # short last step reaches the endpoint
    # The decision is made from the step COUNT, not from a time tolerance:
    # far from the origin an ulp of time can exceed the step itself, and a
    # tolerance-only snap would declare the run over after the first step.
    far = IVP((u, t) -> one.(u), [0.0], (1e16, 1e16 + 8.0))
    sol = solve(far, Euler(h=2.0))
    @test length(sol.t) == 5 && sol.t[end] == 1e16 + 8.0
    @test sol.u[end][1] ≈ 8.0
    # …and a schedule that does NOT divide the span is not declared complete
    # early, however coarse the time grid is out there: h = 5 on a span of 12
    # is 2.4 steps, so the solver must take three, shortening the last step.
    far2 = IVP((u, t) -> one.(u), [0.0], (1e16, 1e16 + 12.0))
    sol2 = solve(far2, Euler(h=5.0))
    @test sol2.u[end][1] ≈ 12.0 && sol2.t[end] == 1e16 + 12.0
    # Adaptive runs also finish exactly on the requested endpoint:
    ad = solve(Logistic(0.3, (0.0, 1.0)), DormandPrince54(h=0.1, εᵣ=1e-6, save_stepsizes=true))
    @test ad.t[end] == 1.0
end

@testset "Explicit chunk schedules" begin
    for dt in (1.6, 500/99, 0.003), h in (0.01, 0.1), method in (RK4, Euler)
        c=chunk_solver(method(h=h), dt)
        p=Logistic([0.1], (0.0, dt))
        a=solve(p, c)
        b=solve(Logistic([0.1], (103.0, 103.0 + dt)), c)
        @test a.u[end]==b.u[end]
        @test a.t[end]==dt
        @test b.t[end]==103.0 + dt
        @test length(a.u)==c.plan.full+(c.plan.tail>0)+1
        @test c.solver.stepsize.h==h
    end
    @test_throws ArgumentError chunk_solver(RK4(h=0.1), 0.0)
    @test_throws ArgumentError solve(Logistic([0.1], (0.0, 2.0)), chunk_solver(RK4(h=0.1), 1.0))
end

@testset "Short final step across solver families" begin
    p = Dahlquist(1.0, (0.0, 0.23))
    for method in (RK4, BackwardEuler, GaussLegendre4, ExponentialEuler)
        solver = method(h=0.1)
        sol = solve(p, solver)
        first = solve(NSDEBase.remake(p; tspan=(0.0, 0.2)), solver)
        tail = solve(
            NSDEBase.remake(p; u0=first.u[end], tspan=(0.2, 0.23)),
            NSDEBase.withstepsize(solver, 0.23-0.2),
        )
        @test sol.t[end] == 0.23
        @test length(sol.t) == 4
        @test sol.u[end] ≈ tail.u[end] rtol=1e-14
        @test solver.stepsize.h == 0.1
    end
    for method in (Euler, DormandPrince54)
        solver = method(h=0.1)
        sol = solve(p, solver; dense=true)
        @test sol.t[end] == 0.23
        @test sol(0.23, solver.tableau) == sol.u[end]
        @test isfinite(sol(0.22, solver.tableau)[1])
    end
    zero_span = NSDEBase.remake(p; tspan=(2.0, 2.0))
    @test solve(zero_span, RK4(h=0.1)).u == [p.u0]
    @test_throws ArgumentError solve(NSDEBase.remake(p; tspan=(1.0, 0.0)), RK4(h=0.1))
end

@testset "Rounding-sized tails at translated endpoints" begin
    span=(15.200000000000001, 16.200000000000003)
    sol=solve(Logistic(0.3, span), RK4(h=1e-3))
    @test length(sol.t)==1001
    @test sol.t[end]==span[2]
    @test all(diff(sol.t) .> 9e-4)
    p=IVP((u, t)->ones(length(u)), [0.0], (6000.0, 6000.6))
    sol=solve(p, RK4(h=1e-4))
    @test length(sol.t)==6001
    @test sol.u[end][1] ≈ p.tspan[2]-p.tspan[1] atol=1e-12
    # Both sides of a dividing duration avoid a near-zero final step.
    for dt in (1.0-1e-12, 1.0+1e-12)
        plan=NSDERungeKutta.StepSchedule(dt, 1e-3)
        @test plan.full+(plan.tail>0)==1000
        @test plan.tail==0 || plan.tail>0.9e-3
    end
    # A genuine, short, nondividing final interval is still integrated.
    sol=solve(IVP((u, t)->ones(length(u)), [0.0], (0.0, 1.0+1e-7)), RK4(h=1e-3))
    @test length(sol.t)==1002
    @test sol.u[end][1] ≈ 1.0+1e-7 atol=1e-14
    for start in (0.0, 15.2, 6000.0)
        p=IVP((u, t)->ones(length(u)), [0.0], (start, start+0.6))
        sol=solve(p, DormandPrince54(h=0.2, εₐ=1e-12, save_stepsizes=true))
        @test sol.t[end]==p.tspan[2]
        @test minimum(sol.stepsizes.accepted)>1e-6
        @test sol.u[end][1] ≈ p.tspan[2]-start atol=2e-12
    end
end

@testset "Cached exponential final-step coefficients" begin
    calls=Ref(0)
    method=ExponentialEuler(h=0.1)
    tab=method.tableau
    wrapped=z->(calls[]+=1; tab.κ(z))
    counted=NSDERungeKutta.ExponentialTableau(tab.name, tab.p, tab.q, tab.s, tab.c, wrapped)
    solver=NSDERungeKutta.EXPRK(counted, 0.1)
    p=Dahlquist(1.0, (0.0, 0.23))
    planned=chunk_solver(solver, 0.23)
    cache=initialize_cache(p, planned)
    sol=initialize_solution(p, planned)
    built=calls[]
    @test built==2
    for start in (0.0, 15.2, 6000.0)
        shifted=remake(p; tspan=(start, start+0.23))
        solve!(cache, sol, shifted, planned)
        @test sol.u[end][1] ≈ exp(0.23) rtol=1e-13
        @test sol.t[end]==shifted.tspan[2]
    end
    @test calls[]==built
end

@testset "Plain exponential solve reuses tail coefficients" begin
    calls = Ref(0)
    tab = ExponentialEuler(h=0.1).tableau
    counted = NSDERungeKutta.ExponentialTableau(
        tab.name, tab.p, tab.q, tab.s, tab.c, z->(calls[]+=1; tab.κ(z))
    )
    solver = NSDERungeKutta.EXPRK(counted, 0.1)
    p = IVP([-2.0 1.0; 0.0 -1.0], [1.0, 1.0], (0.0, 0.23))
    cache = initialize_cache(p, solver)
    solution = initialize_solution(p, solver)
    @test calls[] == 1
    for _ in 1:4
        solve!(cache, solution, p, solver)
        @test solution.u[end] ≈ exp(0.23*p.rhs.L)*p.u0 rtol=1e-13
        @test calls[] == 2
    end
    longer = remake(p; tspan=(0.0, 0.27))
    solve!(cache, solution, longer, solver)
    @test calls[] == 3
    solve!(cache, solution, longer, solver)
    @test calls[] == 3
end
