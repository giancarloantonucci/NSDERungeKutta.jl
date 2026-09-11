using NSDERungeKutta
using LinearAlgebra
using Test
using Aqua
import RecipesBase
import NSDEBase

# RecipesBase declares `is_key_supported` but leaves it for the plotting
# FRONTEND (Plots) to define; headless, recipes with KEYWORD ARGUMENTS hit it
# when deciding whether an extracted kwarg stays in the attribute dict, and
# error without a method. Answering `true` for everything is the standard
# headless-testing shim: it affects only which keys survive in the dict,
# never the recipe logic under test. (Test-only; nothing ships this.)
RecipesBase.is_key_supported(::Symbol) = true

const AP = NSDERungeKutta.AdaptiveParameters
const NP = NSDERungeKutta.NewtonParameters

# ------------------------------------------------------------------ helpers --

# THE abstraction audit: one harness, any solver, any problem, no family
# branch. Solve at h₀, h₀/2, …, take the error at the final grid point against
# the exact solution, and estimate the observed order as the MAX of
#   (a) the least-squares log-log slope, and
#   (b) the finest-grid pairwise ratio (backing off one rung at roundoff).
# The max is deliberate: near a cancellation of the leading error term, or on
# a slow entry to the asymptotic regime, (a) under-reads while (b) is sound;
# at the roundoff floor (b) degrades while (a) holds. A genuinely wrong
# tableau shows its true (lower) order on BOTH estimators — verified against
# deliberately re-broken tableaus during development.
function convergence_order(make_solver, problem, exact; h0, refinements=4)
    hs = [h0 / 2.0^j for j in 0:refinements-1]
    errs = map(hs) do h
        solution = solve(problem, make_solver(h))
        norm(solution.u[end] - exact(solution.t[end])) + 1e-300
    end
    x, y = log.(hs), log.(errs)
    x̄, ȳ = sum(x) / length(x), sum(y) / length(y)
    lsq = sum((x .- x̄) .* (y .- ȳ)) / sum(abs2, x .- x̄)
    j = length(errs)
    if errs[j] < 50 * eps() && length(errs) > 2
        j -= 1 # finest error is machine noise: use the previous rung
    end
    ratio = log(errs[j-1] / errs[j]) / log(hs[j-1] / hs[j])
    return max(lsq, ratio)
end

# h₀ and refinement count by order: high-order methods hit roundoff fast, so
# they get larger steps and fewer halvings to keep every error in the
# resolvable band.
function sweep_params(p)
    p ≤ 2 && return (h0=0.1, refinements=4, below=0.4)
    p ≤ 4 && return (h0=0.2, refinements=4, below=0.4)
    p ≤ 5 && return (h0=0.4, refinements=4, below=0.4)
    return (h0=0.5, refinements=3, below=0.7) # 3-point ladders on p ≥ 6 carry more slack
end

# Under-order is the bug signal (wrong tableau, wrong kernel); measuring a bit
# high on smooth problems is benign. Hence the asymmetric window.
function assert_order(name, make_solver, problem, exact, p)
    prm = sweep_params(p)
    slope = convergence_order(make_solver, problem, exact; h0=prm.h0, refinements=prm.refinements)
    @test p - prm.below < slope < p + 2.5
    if !(p - prm.below < slope < p + 2.5)
        @info "order mismatch" name p slope
    end
end

function step_allocs(problem, solver)
    solution = RungeKuttaSolution(problem, solver)
    cache = NSDERungeKutta.RungeKuttaCache(problem, solver)
    NSDERungeKutta.step!(cache, solution, problem, solver) # warm-up/compile
    return @allocated NSDERungeKutta.step!(cache, solution, problem, solver)
end

# ------------------------------------------------------------------ problems --

dahlquist = Dahlquist(0.5, (0.0, 1.0))                    # LinearRightHandSide
dahlquist_exact = t -> [0.5 * exp(t)]
logistic = Logistic(0.3, (0.0, 1.0))                      # NonlinearRightHandSide
logistic_exact = t -> [0.3 * exp(t) / (1.0 + 0.3 * (exp(t) - 1.0))]
split_lin = IVP(hcat(-0.5), (u, t) -> 1.5 .* u, [0.3], (0.0, 1.0))            # LRHS stiff part
split_nlin = IVP(SRHS((u, t) -> -0.5 .* u, (u, t) -> 1.5 .* u), [0.3], (0.0, 1.0)) # NRHS stiff part
split_exact = t -> [0.3 * exp(t)]

# ------------------------------------------------------------------ the zoo --

erk_fixed = [Euler, Heun2, Midpoint, Ralston2, Heun3, RungeKutta3, Ralston3,
             SSPRK3, Ralston4, RungeKutta4, Rule38, Butcher5, KuttaNystrom5,
             Butcher6, Butcher7]
erk_embedded = [HeunEuler, BogackiShampine, Fehlberg45, DormandPrince54,
                Verner65, Fehlberg78]
dirk_zoo = [BackwardEuler, ImplicitMidpoint, SDIRK2, LobattoIII2, CrankNicolson,
            SDIRK3, RadauI3, RadauII3, SDIRK4, LobattoIII4]
irk_zoo = [LobattoIIIC2, RadauIA3, RadauIIA3, GaussLegendre4, LobattoIIIA4,
           LobattoIIIB4, LobattoIIIC4, RadauI5, RadauIA5, RadauII5, RadauIIA5,
           GaussLegendre6]
ierk_zoo = [IMEXEuler, IMEXSSP2_222, IMEXSSP2_322, IMEXSSP2_332, IMEXSSP3_332]

@testset "NSDERungeKutta" begin

@testset "Aqua" begin
    # Same reasoning as NSDEBase: recipetype(::Val{…}) is the established
    # user-plot idiom; whitelist that one function, keep the rest strict.
    Aqua.test_all(NSDERungeKutta; piracies=(; treat_as_own=[NSDERungeKutta.RecipesBase.recipetype]))
end

@testset "Butcher tableau sanity" begin
    tableaus = Tuple{String, Any}[]
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

@testset "adaptive stepping (explicit-only in v1)" begin
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
        hs = recording.stepsize.hs
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
        # boundary made visible. (From the thesis suite's
        # AdaptiveStepsizeDemo.jl.)
        uexact(t) = [2exp(-t) + sin(t), 2exp(-t) + cos(t)]
        f₁(u, t) = [-2.0 1.0; 1.0 -2.0] * u + [2sin(t), 2(cos(t) - sin(t))]        # {-1, -3}
        f₂(u, t) = [-2.0 1.0; 998.0 -999.0] * u + [2sin(t), 999(cos(t) - sin(t))]  # {-1, -1000}
        function adaptive_run(f)
            solver = Fehlberg45(h=0.1, εᵣ=1e-2, save_stepsizes=true)
            sol = solve(IVP(RHS(f), [2.0, 3.0], (0.0, 10.0)), solver)
            nsteps = numtimesteps(sol) - 1
            nrejections = sum(length, solver.stepsize.hs.rejected)
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

@testset "dense output and interpolation" begin
    @testset "stage history aligns with accepted steps" begin
        problem = Logistic(0.3, (0.0, 1.0))
        solver = BogackiShampine(h=0.05, εᵣ=1e-6)
        solution = solve(problem, solver; dense=true)
        @test length(solution.k) == numtimesteps(solution) - 1 # the drift fix
        for tₚ in (0.13, 0.5, 0.87)
            @test norm(solution(tₚ, solver.tableau) - logistic_exact(tₚ)) < 1e-4
        end
    end

    @testset "spline fallbacks" begin
        solution = solve(dahlquist, RungeKutta4(h=0.1))
        exact = dahlquist_exact(0.55)
        err_linear = norm(solution(0.55) - exact)
        err_hermite = norm(solution(0.55, (u, t) -> u) - exact) # u′ = u
        @test err_linear < 5e-3
        @test err_hermite < 1e-5
        @test err_hermite < err_linear
        @test solution(-10.0) == solution.u[1]   # clamped below
        @test solution(99.0) == solution.u[end]  # clamped above
        @test solution(solution.t[3]) == solution.u[3] # exact at grid nodes (pins the segment lookup)
    end
end

@testset "allocation-free ERK step" begin
    L = [0.0 1.0; -4.0 -0.4]
    g!(dg, t) = (dg[1] = sin(t); dg[2] = cos(t); dg)
    f!(du, u, t) = (@. du = -u; du)
    problems = (
        IVP(f!, ones(2), (0.0, 1.0)),                                     # NonlinearRHS
        IVP(LRHS(L, g!), ones(2), (0.0, 1.0)),                            # LRHS with forcing
        IVP(SRHS(LRHS(L, g!), RHS((du, u, t) -> (@. du = sin(u); du))), ones(2), (0.0, 1.0)), # split
    )
    for problem in problems
        @test step_allocs(problem, RungeKutta4(h=0.5)) == 0
    end
end

@testset "solution container" begin
    solution = solve(dahlquist, RungeKutta4(h=0.25))
    @test numvariables(solution) == 1
    @test numtimesteps(solution) == length(solution.t)
    @test length(extract(solution, 1)) == numtimesteps(solution)

    # setindex! with a solution slice (the fixed `@↓ … ← …` path)
    solution[2] = solution[1]
    @test solution.u[2] == solution.u[1]
    @test solution.t[2] == solution.t[1]
    # tuple form
    solution[3] = ([9.9], 9.9)
    @test solution.u[3] == [9.9] && solution.t[3] == 9.9
end

@testset "smoke: solve runs and returns the right type" begin
    for solver in (RungeKutta4(h=0.1), BackwardEuler(h=0.1), GaussLegendre4(h=0.1),
                   IMEXEuler(h=0.1), DormandPrince54(h=0.1))
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

@testset "RecipesBase recipes (headless)" begin
    # `apply_recipe` exercises the Plots-side recipe surface with no Plots
    # dependency — previously the one exported surface the suite never
    sol = solve(Lorenz([2.0, 3.0, -14.0], (0.0, 1.0)), RK4(h=1e-2))
    for attrs in (Dict{Symbol,Any}(),
                  Dict{Symbol,Any}(:variables => 2),
                  Dict{Symbol,Any}(:iscomplex => false, :skip => 10))
        @test RecipesBase.apply_recipe(attrs, sol) isa Vector{RecipesBase.RecipeData}
    end
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(), NSDEBase._PhasePlot(sol)) isa Vector{RecipesBase.RecipeData}
    zsol = solve(Dahlquist([1.0 + 0.0im], (0.0, 1.0); λ=-0.5 + 1.0im), RK4(h=1e-2))
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:iscomplex => true), zsol) isa Vector{RecipesBase.RecipeData}
    # The step-size record recipes, including the trailing-bin guard:
    recording = DormandPrince54(h=0.1, εᵣ=1e-4, save_stepsizes=true)
    rsol = solve(Logistic(0.3, (0.0, 1.0)), recording)
    hs = recording.stepsize.hs
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(), hs) isa Vector{RecipesBase.RecipeData}
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(), rsol.t[1:end-1], hs) isa Vector{RecipesBase.RecipeData}
end


@testset "a dividing step ends exactly on tN" begin
    # 0.1 is not representable, so three steps of 0.1 from 0 need not land on
    # 0.3 bit for bit; a landing within rounding is snapped onto tN so the
    # solve neither takes a spurious fourth step nor stops an ulp past the end.
    for (tspan, h, nsteps) in (((0.0, 0.3), 0.1, 3), ((0.0, 1.0), 1e-3, 1000), ((0.25, 0.5), 1e-3, 250), ((0.0, 0.7), 0.1, 7))
        sol = solve(Logistic(0.3, tspan), RK4(h=h))
        @test sol.t[end] == tspan[2]
        @test length(sol.t) == nsteps + 1
    end
    over = solve(Logistic(0.3, (0.0, 0.3)), RK4(h=0.2))     # 0.2 does not divide 0.3
    @test over.t[end] ≈ 0.4 && length(over.t) == 3          # honest overshoot, untouched
    # The decision is made from the step COUNT, not from a time tolerance:
    # far from the origin an ulp of time can exceed the step itself, and a
    # tolerance-only snap would declare the run over after the first step.
    far = IVP((u, t) -> one.(u), [0.0], (1e16, 1e16 + 8.0))
    sol = solve(far, Euler(h=2.0))
    @test length(sol.t) == 5 && sol.t[end] == 1e16 + 8.0
    @test sol.u[end][1] ≈ 8.0
    # …and a schedule that does NOT divide the span is not declared complete
    # early, however coarse the time grid is out there: h = 5 on a span of 12
    # is 2.4 steps, so the solver must take three and overshoot as documented.
    far2 = IVP((u, t) -> one.(u), [0.0], (1e16, 1e16 + 12.0))
    @test NSDERungeKutta.fixedstep_count(Euler(h=5.0), 1e16, 1e16 + 12.0) == 0
    sol2 = solve(far2, Euler(h=5.0))
    @test sol2.u[end][1] ≈ 15.0 && sol2.t[end] ≥ 1e16 + 12.0
    @test NSDERungeKutta.fixedstep_count(Euler(h=1e-3), 0.0, 1.0) == 1000
    @test NSDERungeKutta.fixedstep_count(Euler(h=1e-3), 0.25, 0.5) == 250
    @test NSDERungeKutta.fixedstep_count(Euler(h=0.7), 0.0, 0.3) == 0  # fewer than one step
    @test NSDERungeKutta.fixedstep_count(RK4(h=0.1), 0.0, 0.3) == 3
    @test NSDERungeKutta.fixedstep_count(RK4(h=0.2), 0.0, 0.3) == 0
    @test NSDERungeKutta.fixedstep_count(DormandPrince54(h=0.1), 0.0, 0.3) == 0 # adaptive: no schedule
    # Adaptive runs are not touched (their end handling is the controller's):
    ad = solve(Logistic(0.3, (0.0, 1.0)), DormandPrince54(h=0.1, εᵣ=1e-6, save_stepsizes=true))
    @test ad.t[end] ≥ 1.0
end

@testset "dense slicing: stages live on intervals, not nodes" begin
    # A dense solution with N nodes stores N − 1 stage sets. The indexers used
    # to read `k[i]` at node i, so `sol[end]` threw a BoundsError.
    dense = solve(dahlquist, RungeKutta4(h=0.5); dense=true)   # t = 0, 0.5, 1
    N = length(dense)
    @test length(dense.k) == N - 1
    last = dense[end]                                          # used to throw
    @test last.u == [dense.u[end]] && last.t == [dense.t[end]] && last.k === nothing
    @test dense[1].k === nothing                               # a node owns no interval
    slice = dense[1:2]                                         # two nodes, one interval
    @test slice.u == dense.u[1:2] && length(slice.k) == 1 && slice.k[1] == dense.k[1]
    whole = dense[1:N]
    @test length(whole.k) == N - 1                             # `end` inside a range is fine too
    @test_throws ArgumentError dense[[1, N]]                   # non-adjacent: no stages between them
    sparse = solve(dahlquist, RungeKutta4(h=0.5))
    @test sparse[[1, N]].k === nothing                         # non-dense: any selection is fine
    # Writes into a DENSE solution: partial writes are refused BEFORE any
    # mutation (the adjoining intervals' stages could no longer match), a
    # whole-solution write with matching stages goes through.
    before = (copy.(dense.u), copy(dense.t), deepcopy(dense.k))
    @test_throws ArgumentError dense[2] = dense[1]
    @test_throws ArgumentError dense[2] = ([9.9], 9.9)
    @test_throws ArgumentError dense[1:2] = slice
    @test_throws ArgumentError dense[1:N] = sparse                # nodes only: no stages to go with them
    @test_throws DimensionMismatch dense[1:2] = dense[1:N]        # sizes disagree
    @test dense.u == before[1] && dense.t == before[2] && dense.k == before[3] # untouched by every refusal
    other = solve(Dahlquist(0.5, (0.0, 1.0); λ=-2.0), RungeKutta4(h=0.5); dense=true)
    dense[1:N] = other
    @test dense.u == other.u && dense.t == other.t && dense.k == other.k
    # Non-dense targets accept node writes as before:
    sparse[2] = sparse[1]
    @test sparse.u[2] == sparse.u[1]
    u1, u3 = copy(sparse.u[1]), copy(sparse.u[3])
    sparse[[1, 3]] = sparse[[3, 1]]                               # non-contiguous is fine without stages
    @test sparse.u[1] == u3 && sparse.u[3] == u1
end

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
    err = try; solve(blowup, BackwardEuler(h=1.0, Mₙ=2)); nothing; catch e; e; end
    @test err isa NewtonFailure && err.t == 0.0 && err.stage == 1 && err.updates == 2
    @test !(err.residual ≤ err.tolerance)  # holds for Inf/NaN residuals too
    @test occursin("did not converge", sprint(showerror, err))
    # Away from the singularity the same problem steps fine (u = 1/(1 − t)):
    @test solve(IVP((u, t) -> u .^ 2, [1.0], (0.0, 0.5)), BackwardEuler(h=0.005)).u[end][1] ≈ 2.0 rtol = 2e-2
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
        for solver in (RadauIIA3(h=1.0, εₐ=0.0, εᵣ=1e-8), BackwardEuler(h=1.0, εₐ=0.0, εᵣ=1e-8), SDIRK2(h=1.0, εₐ=0.0, εᵣ=1e-8))
            @test solve(tiny, solver).u[end][1] ≈ c rtol = 1e-8
        end
    end
    # Accepted stages satisfy the stage equation to the stated tolerance: for
    # backward Euler, ‖f(u₁) − (u₁ − u₀)/h‖ ≤ εₐ + εᵣ·scale at every step.
    εᵣ, εₐ = 1e-8, 1e-12
    sol = solve(logistic, BackwardEuler(h=0.1, εᵣ=εᵣ, εₐ=εₐ, Mₙ=20))
    for n = 1:length(sol) - 1
        k = (sol.u[n+1] - sol.u[n]) / 0.1
        f = logistic.rhs(sol.u[n+1], sol.t[n+1])
        @test norm(f - k) ≤ εₐ + εᵣ * max(norm(k), norm(f))
    end
    @test NP(εᵣ=1e-3, εₐ=1e-9, Mₙ=4).εₐ == 1e-9
    @test NP(1e-3, 10).εₐ == 1e-12                                     # legacy positional form
    @test_throws ArgumentError NP(εᵣ=-1.0)
    @test_throws ArgumentError NP(Mₙ=0)
end

end # outer testset

include("exprk.jl")
include("oracle.jl") # mathematical oracle for the exponential schemes; needs Printf (test extra)
