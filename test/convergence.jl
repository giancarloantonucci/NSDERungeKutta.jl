@testset "convergence order: fixed-step ERK" begin
    for F in erk_fixed
        p = F(h=0.1).tableau.p
        for (problem, exact) in ((dahlquist, dahlquist_exact), (logistic, logistic_exact))
            assert_order(string(F), h -> F(h=h), problem, exact, p)
        end
    end
end

@testset "convergence order: embedded tableaus run fixed-step" begin
    # Strip the adaptive controller and march the b-weights at fixed h: the
    # slope must match the tableau's claimed p (catches b/p mislabels).
    for F in erk_embedded
        tableau = F(h=0.1).tableau
        assert_order(string(F), h -> ERK(tableau, h), dahlquist, dahlquist_exact, tableau.p)
    end
end

@testset "convergence order: DIRK (Newton and direct-linear branches)" begin
    for F in dirk_zoo
        p = F(h=0.1).tableau.p
        make = h -> F(h=h, εᵣ=1e-12, Mₙ=30)
        assert_order(string(F), make, dahlquist, dahlquist_exact, p)  # LRHS: direct solve
        assert_order(string(F), make, logistic, logistic_exact, p)    # NRHS: simplified Newton
    end
end

@testset "convergence order: IRK (coupled solve)" begin
    for F in irk_zoo
        p = F(h=0.1).tableau.p
        make = h -> F(h=h, εᵣ=1e-12, Mₙ=30)
        assert_order(string(F), make, dahlquist, dahlquist_exact, p)
        assert_order(string(F), make, logistic, logistic_exact, p)
    end
end

@testset "convergence order: IMEX on split problems" begin
    for F in ierk_zoo
        p = F(h=0.1).implicitableau.p
        make = h -> F(h=h, εᵣ=1e-12, Mₙ=30)
        assert_order(string(F), make, split_lin, split_exact, p)   # linear stiff part
        assert_order(string(F), make, split_nlin, split_exact, p)  # nonlinear stiff part
    end
end

@testset "smoke: solve runs and returns the right type" begin
    for solver in (
        RungeKutta4(h=0.1),
        BackwardEuler(h=0.1),
        GaussLegendre4(h=0.1),
        IMEXEuler(h=0.1),
        DormandPrince54(h=0.1),
    )
        problem = solver isa ImplicitExplicitRungeKuttaSolver ? split_lin : logistic
        @test solve(problem, solver) isa RungeKuttaSolution
    end
end

@testset "stiff problem: implicit stable where explicit blows up" begin
    # The implicit families exist for stiffness, so prove it: hλ = -100 sits
    # far outside every explicit stability region and comfortably inside the
    # A-stable ones. The exact solution decays to ~0 instantly.
    stiff = Dahlquist(1.0, (0.0, 1.0); λ=-1e5)
    h = 1e-3
    explicit = solve(stiff, RungeKutta4(h=h))
    u_exp = norm(explicit.u[end])
    @test !isfinite(u_exp) || u_exp > 1e10                     # RK4 detonates
    for solver in (BackwardEuler(h=h), CrankNicolson(h=h), RadauIIA3(h=h, εᵣ=1e-12, Mₙ=30))
        implicit = solve(stiff, solver)
        @test all(isfinite, implicit.u[end])
        @test norm(implicit.u[end]) ≤ 1.0                      # |R(hλ)| ≤ 1: no growth, ever
    end
    # And the L-stable one actually kills the transient rather than ringing:
    @test norm(solve(stiff, BackwardEuler(h=h)).u[end]) < 1e-10
end
