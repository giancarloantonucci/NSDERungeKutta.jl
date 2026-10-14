@testset "Newton: no silent use of an unconverged stage" begin
    # u′ = u², u₀ = 1, backward Euler with h = 1: the stage equation
    # y = 1 + y² has no real root. The solver used to return a finite garbage
    # step after exhausting Mₙ, with no signal at all; a later draft compared
    # `Inf ≤ Inf` and returned -Inf. Neither is acceptable.
    blowup = IVP((u, t) -> u .^ 2, [1.0], (0.0, 1.0))
    @test_throws NewtonFailure solve(blowup, BackwardEuler(h=1.0, Mₙ=2))
    @test_throws NewtonFailure solve(blowup, BackwardEuler(h=1.0, Mₙ=50)) # diverges to ±Inf: still throws
    @test_throws NewtonFailure solve(blowup, BackwardEuler(h=1.0, Mₙ=1000))
    # A tiny UPDATE is not convergence. Here the first simplified-Newton
    # update is 1e-20 while the stage residual is about 1e20: an
    # increment-based test with any absolute tolerance accepts it.
    spiky = IVP((u, t) -> @.(1 - 1e20 * u + 1e60 * u^2), [0.0], (0.0, 1.0))
    @test_throws NewtonFailure solve(spiky, BackwardEuler(h=1.0))
    @test_throws NewtonFailure solve(spiky, BackwardEuler(h=1.0, εₐ=1e-6, εᵣ=1e-2))
    # Every implicit family has the same exit. One update from a zero start
    # cannot satisfy a tight residual bound on a nonlinear stage equation:
    tight = (Mₙ=1, εᵣ=1e-14, εₐ=0.0)
    @test_throws NewtonFailure solve(logistic, RadauIIA3(h=0.5; tight...))   # IRK
    split = IVP((u, t) -> -0.3 .* u .^ 2, (u, t) -> 0.3 .* u, [0.5], (0.0, 1.0)) # nonlinear stiff part
    @test_throws NewtonFailure solve(split, IMEXEuler(h=0.5; tight...))      # IERK
    @test_throws NewtonFailure solve(logistic, SDIRK2(h=0.5; tight...))       # DIRK, later stage
    @test solve(logistic, SDIRK2(h=0.5)).u[end] ≈ logistic_exact(1.0) atol = 5e-2 # defaults still converge
    # The exception says where and how badly:
    err = try
        solve(blowup, BackwardEuler(h=1.0, Mₙ=2))
        nothing
    catch e
        e
    end
    @test err isa NewtonFailure && err.t == 0.0 && err.stage == 1 && err.updates == 2
    @test !(err.residual ≤ err.tolerance)  # holds for Inf/NaN residuals too
    @test occursin("did not converge", sprint(showerror, err))
    # Away from the singularity the same problem steps fine (u = 1/(1 − t)):
    @test solve(IVP((u, t) -> u .^ 2, [1.0], (0.0, 0.5)), BackwardEuler(h=0.005)).u[end][1] ≈ 2.0 rtol =
        2e-2
    # Non-finite tolerances are refused, so `Inf ≤ Inf` can never be the test:
    @test_throws ArgumentError NP(εᵣ=Inf)
    @test_throws ArgumentError NP(εₐ=NaN)
    @test_throws ArgumentError BackwardEuler(h=0.1, εᵣ=Inf)
end

@testset "Newton: residual acceptance, not increment size" begin
    # At an exact equilibrium the stage residual is identically zero and the
    # stage is accepted with no update at all (a purely relative increment
    # test could never pass there, and with a failure exit would throw):
    rest = IVP((u, t) -> u .* (1 .- u), [1.0], (0.0, 1.0))
    for solver in (BackwardEuler(h=0.5, Mₙ=3), RadauIIA3(h=0.5, Mₙ=3), SDIRK3(h=0.5, Mₙ=3))
        solution = solve(rest, solver)
        @test all(u -> u == [1.0], solution.u)
    end
    # For a split problem the equilibrium must be one of BOTH parts: the
    # splitting does not preserve a state at which the two vector fields
    # merely cancel each other (u′ = −u² + u at u = 1 is stepped to ≈ 1.098
    # by IMEX Euler, correctly for that method).
    restsplit = IVP((u, t) -> 0.3 .* (1 .- u), (u, t) -> (1 .- u) .^ 2, [1.0], (0.0, 1.0))
    @test all(u -> u == [1.0], solve(restsplit, IMEXEuler(h=0.5, Mₙ=3)).u)
    # Residual norms must not under- or overflow for representable scales:
    # the coupled IRK residual used to be summed as squares, so u′ = 1e-200
    # was "accepted" at the zero initial stages and u′ = 1e200 threw.
    for c in (1e-200, 1e-100, 1.0, 1e100, 1e200)
        tiny = IVP((u, t) -> fill(c, length(u)), [0.0], (0.0, 1.0))
        for solver in (
            RadauIIA3(h=1.0, εₐ=0.0, εᵣ=1e-8),
            BackwardEuler(h=1.0, εₐ=0.0, εᵣ=1e-8),
            SDIRK2(h=1.0, εₐ=0.0, εᵣ=1e-8),
        )
            @test solve(tiny, solver).u[end][1] ≈ c rtol = 1e-8
        end
    end
    # Accepted stages satisfy the stage equation to the stated tolerance: for
    # backward Euler, ‖f(u₁) − (u₁ − u₀)/h‖ ≤ εₐ + εᵣ·scale at every step.
    εᵣ, εₐ = 1e-8, 1e-12
    sol = solve(logistic, BackwardEuler(h=0.1, εᵣ=εᵣ, εₐ=εₐ, Mₙ=20))
    for n in 1:(length(sol) - 1)
        k = (sol.u[n + 1] - sol.u[n]) / 0.1
        f = logistic.rhs(sol.u[n + 1], sol.t[n + 1])
        @test norm(f - k) ≤ εₐ + εᵣ * max(norm(k), norm(f))
    end
    @test NP(εᵣ=1e-3, εₐ=1e-9, Mₙ=4).εₐ == 1e-9
    @test NP(1e-3, 10).εₐ == 1e-12                                     # legacy positional form
    @test_throws ArgumentError NP(εᵣ=-1.0)
    @test_throws ArgumentError NP(Mₙ=0)
end
