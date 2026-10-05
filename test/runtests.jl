using NSDEMovingWindow
using NSDERungeKutta
using LinearAlgebra
using Test
using Aqua
import RecipesBase
import NSDEBase

# Truth for accuracy assertions: the serial fine solve over the whole span.
const problem = Logistic(0.3, (0.0, 3.0))
const finesolver = RK4(h=1e-3)
const finetruth = solve(problem, finesolver)

good_parareal() = Parareal(finesolver, RungeKutta4(h=5e-2);
                           parameters=PararealParameters(N=4, K=4), tolerance=Tolerance(ϵ=1e-9))
# A setup built to FAIL windows: hopeless coarse solver, tiny budget, tight ϵ.
# ϵ sits above what a shrunken/loosened window reaches after a restart or two,
# so the adaptive machinery both FIRES and RECOVERS — a tolerance far below
# reach would send Stretch/Leap into a shrink spiral (Δτ collapsing, window
# count exploding), which is the algorithm as designed, not a test's place.
bad_parareal(; ψ=ψ₁, weights=Weights()) =
    Parareal(finesolver, Euler(h=0.5);
             parameters=PararealParameters(N=4, K=2), tolerance=Tolerance(ϵ=1e-5, ψ=ψ, weights=weights))

windowstarts(sol) = [NSDEMovingWindow.windowspan(sol[m])[1] for m = 1:length(sol)]
function windowlengths(sol)
    return [let (ta, tb) = NSDEMovingWindow.windowspan(sol[m]); tb - ta; end for m = 1:length(sol)]
end

@testset "NSDEMovingWindow" begin

@testset "Aqua" begin
    Aqua.test_all(NSDEMovingWindow)
end

@testset "construction validation" begin
    parareal = good_parareal()
    @test_throws ArgumentError MoWi(parareal; τ=1.0, Δτ=2.0)   # Δτ > τ
    @test_throws ArgumentError MoWi(parareal; τ=1.0, Δτ=0.0)   # no forward shift
    @test_throws ArgumentError MoWi(parareal; τ=0.0)           # no window
    # Zoom over the unweighted ψ₁ would be a silent no-op: refuse loudly.
    @test_throws ArgumentError MoWi(parareal; τ=1.0, adaptive=ZoomParameters())
    weighted = bad_parareal(ψ=ψ₂, weights=Weights(w=2.0))
    @test MoWi(weighted; τ=1.0, adaptive=ZoomParameters()) isa MoWi
    # Zoom + updatew: `update!`'s max(w, Λ̂) floor undoes every δw⁺ tightening
    # on the next sweep — refuse the pairing loudly at construction.
    remeasured = bad_parareal(ψ=ψ₂, weights=Weights(updatew=true))
    @test_throws ArgumentError MoWi(remeasured; τ=1.0, adaptive=ZoomParameters())
    @test MoWi(remeasured; τ=1.0) isa MoWi                          # fine without Zoom
    @test MoWi(remeasured; τ=1.0, adaptive=LeapParameters()) isa MoWi # Leap never touches w
    # The solver is deep-copied: adaptive runs must not mutate the user's.
    mowi = MoWi(parareal; τ=1.0)
    @test mowi.parallelsolver !== parareal
end

@testset "flatten" begin
    mowi = MoWi(good_parareal(); τ=1.0, Δτ=0.5)
    solution = solve(problem, mowi)
    u, t = flatten(solution)
    @test length(u) == length(t)
    @test length(u) == sum(length(NSDETimeParallel.flatten(solution[m])[1]) for m = 1:length(solution))
    @test !issorted(t) || length(solution) == 1 # window seams rewind time by design
end

@testset "fixed MoWi: end-to-end against the fine solve" begin
    mowi = MoWi(good_parareal(); τ=1.0, Δτ=0.5)
    solution = solve(problem, mowi)
    @test length(solution) == 5                                # (3−1)/0.5 + 1
    @test all(==(0), solution.restarts)
    @test windowstarts(solution) ≈ 0.0:0.5:2.0                 # boundaries advance by Δτ
    @test NSDEMovingWindow.windowspan(solution[end])[2] ≥ 3.0  # full coverage
    for t in 0.0:0.3:3.0
        @test norm(solution(t) - finetruth(t)) < 1e-6
    end
    # Overlap consistency: neighbouring windows agree on their shared stretch.
    for t in 0.5:0.1:1.0
        @test norm(solution[2](t) - solution[1](t)) < 1e-6
    end
end

@testset "implicit coarse solver through MoWi: success and failure" begin
    # Success: backward Euler as the coarse propagator of every window.
    implicit = Parareal(finesolver, BackwardEuler(h=5e-2);
                        parameters=PararealParameters(N=4, K=4), tolerance=Tolerance(ϵ=1e-9))
    solution = solve(problem, MoWi(implicit; τ=1.0, Δτ=0.5))
    @test length(solution) == 5
    @test all(==(0), solution.restarts)
    for t in 0.0:0.3:3.0
        @test norm(solution(t) - finetruth(t)) < 1e-6
    end
    # Failure: a coarse Newton failure inside a window is not something the
    # window loop may swallow or count as a restart; it reaches the caller.
    blow = IVP((u, t) -> u .^ 2, [1.0], (0.0, 0.9))
    failing = Parareal(RK4(h=1e-3), BackwardEuler(h=0.3);
                       parameters=PararealParameters(N=3, K=3), tolerance=Tolerance(ϵ=1e-9))
    @test_throws NewtonFailure solve(blow, MoWi(failing; τ=0.9, Δτ=0.3))
end

@testset "complex states through MoWi" begin
    # KS pushes complex spectral states through the whole MoWi stack in the
    # thesis; CI pins the support here on a complex Dahlquist.
    λ = -0.3 + 1.5im
    zproblem = Dahlquist([1.0 + 0.0im], (0.0, 3.0); λ)
    zfine = solve(zproblem, finesolver)
    zparareal = Parareal(finesolver, RungeKutta4(h=5e-2);
                         parameters=PararealParameters(N=4, K=4), tolerance=Tolerance(ϵ=1e-10))
    zsolution = solve(zproblem, MoWi(zparareal; τ=1.0, Δτ=0.5))
    u, _ = flatten(zsolution)
    @test eltype(first(u)) == ComplexF64                 # nothing silently realified
    for t in 0.0:0.4:3.0
        @test norm(zsolution(t) - zfine(t)) < 1e-7       # window hand-offs preserve ℂ accuracy
    end
end

@testset "strategies on a chaotic problem (smoke)" begin
    # Stretch/Leap/Zoom are pinned mechanically on Logistic with bad_parareal
    # above; this runs all three END-TO-END on Lorenz — the regime the thesis
    # actually deploys them in — at a CI-sized span.
    lproblem = Lorenz([2.0, 3.0, -14.0], (0.0, 6.0))
    ltruth = solve(lproblem, finesolver)
    lparareal(ψ, weights) = Parareal(finesolver, RungeKutta4(h=5e-2);
                                     parameters=PararealParameters(N=4, K=4),
                                     tolerance=Tolerance(ϵ=1e-10, ψ=ψ, weights=weights))
    runs = (
        ("Stretch", MoWi(lparareal(ψ₁, Weights()); τ=2.0, Δτ=1.0, adaptive=StretchParameters())),
        ("Leap",    MoWi(lparareal(ψ₁, Weights()); τ=2.0, Δτ=1.0, adaptive=LeapParameters())),
        # Zoom needs ψ₂ with a FROZEN base weight (the constructor refuses
        # updatew); w = exp(0.9056), the Lorenz Lyapunov rate the thesis uses.
        ("Zoom",    MoWi(lparareal(ψ₂, Weights(w=exp(0.9056))); τ=2.0, Δτ=1.0, adaptive=ZoomParameters())),
    )
    for (name, mowi) in runs
        solution = solve(lproblem, mowi)
        @test NSDEMovingWindow.windowspan(solution[end])[2] ≥ 6.0  # full coverage
        @test all(≤(mowi.adaptive.R), solution.restarts)           # budget respected
        for t in 0.0:1.5:6.0
            # Pointwise tracking on a chaotic system, honestly bounded. The
            # per-window residual (∼ϵ scale) is amplified by LOCAL expansion
            # rates — Lorenz's ‖J‖ bursts far exceed the mean rate Λ ≈ 0.9,
            # so a mean-rate estimate undershoots by orders of magnitude
            # (measured: up to ∼4e-2 absolute by t = 6 at ϵ = 1e-10). The
            # scaled bound below asserts "still on the TRUE trajectory to
            # ∼5%" — which uncontrolled divergence, O(attractor diameter),
            # fails by more than an order of magnitude — with ×10-plus
            # headroom over the measured worst for cross-platform rounding.
            @test norm(solution(t) - ltruth(t)) < 5e-2 * (1 + norm(ltruth(t)))
        end
    end
end

@testset "kwargs reach the subroutine" begin
    mowi = MoWi(good_parareal(); τ=1.5, Δτ=1.5)
    solution = solve(problem, mowi; saveiterates=true)
    @test solution[1].iterates !== nothing                     # the old dispatcher dropped this
end

@testset "Stretch (AMoWi-1): restarts fire, τ shrinks, R caps" begin
    mowi = MoWi(bad_parareal(); τ=1.0, Δτ=0.5, adaptive=StretchParameters(R=3))
    solution = solve(problem, mowi)
    @test any(>(0), solution.restarts)
    @test all(≤(3), solution.restarts)                         # R caps every window
    @test length(unique(round.(windowlengths(solution); digits=10))) > 1 # τ demonstrably moved
    @test NSDEMovingWindow.windowspan(solution[end])[2] ≥ 3.0  # still covers the span
end

@testset "Leap (AMoWi-2): hooks, exemption, cap" begin
    # Unit level (a forced-failure end-to-end can shrink-spiral by design, so
    # the restart mechanics are pinned here deterministically). Floors are set
    # NON-BINDING (0.0) on purpose: this testset pins the hop/leap arithmetic;
    # the floor behaviour has its own testset below.
    state = NSDEMovingWindow.MoWiState(1.0, 0.5, 0.0, 0.0)
    NSDEMovingWindow.adapt_failure!(LeapParameters(), state, nothing)
    @test state.Δτ ≈ 0.125                                     # hop: ×δΔτ⁻
    NSDEMovingWindow.adapt_success!(LeapParameters(δΔτ⁺=100.0), state, nothing)
    @test state.Δτ ≈ 1.0                                       # leap, capped at τ
    @test NSDEMovingWindow.restart_start(LeapParameters(), 7.0, 2.0, 0.25) ≈ 2.25 # retry the SHIFT
    @test NSDEMovingWindow.restart_start(StretchParameters(), 7.0, 2.0, 0.25) ≈ 7.0 # same window start
    @test !NSDEMovingWindow.allowsrestart(LeapParameters(), 1)  # window 1 exempt (thesis §5.2)
    @test NSDEMovingWindow.allowsrestart(LeapParameters(), 2)
    @test NSDEMovingWindow.allowsrestart(StretchParameters(), 1)
    @test !NSDEMovingWindow.allowsrestart(nothing, 2)           # fixed MoWi never restarts
    # Benign setup with aggressive growth: the leap saturates at τ (thesis §5.2/KS).
    mowi = MoWi(good_parareal(); τ=1.0, Δτ=0.25, adaptive=LeapParameters(δΔτ⁺=4.0))
    solution = solve(problem, mowi)
    starts = windowstarts(solution)
    @test length(starts) ≥ 3
    @test starts[3] - starts[2] ≈ 1.0 atol = 1e-12             # capped shift ≡ τ
end

@testset "failure hooks floor their state (termination guards)" begin
    # Without floors, Δτ (Leap) and τ (Stretch) decay geometrically across
    # consecutively failing windows: the advance per window tends to zero and
    # mowi! never terminates. The floors are guards, not adaptive policy.
    parareal = good_parareal()
    N = parareal.parameters.N
    mkstate(parareal) = NSDEMovingWindow.MoWiState(1.0, 0.5, 1.0 / N,
        max(N * parareal.finesolver.stepsize.h, 1.0 / N))
    st = mkstate(parareal)
    for _ = 1:100
        NSDEMovingWindow.adapt_failure!(LeapParameters(δΔτ⁻=0.1), st, parareal)
    end
    @test st.Δτ ≥ st.Δτmin             # floored at one ORIGINAL chunk length
    parareal2 = good_parareal()        # Stretch mutates the coarse step: fresh solver
    st = mkstate(parareal2)
    for _ = 1:100
        NSDEMovingWindow.adapt_failure!(StretchParameters(δτ⁻=0.1), st, parareal2)
    end
    @test st.τ ≥ st.τmin               # window floor holds
    @test st.Δτ ≥ st.Δτmin             # ABSOLUTE shift floor: τ₀/N, not the
    # shrunken τ/N — otherwise the advance decays to one fine step per window
end

@testset "Leap: terminates under persistent failure" begin
    parareal = bad_parareal()
    N = parareal.parameters.N
    tN = problem.tspan[2]
    mowi = MoWi(parareal; τ=1.0, Δτ=0.5, adaptive=LeapParameters(δΔτ⁻=0.1, R=2))
    solution = solve(problem, mowi)    # failures fire; the run must still finish
    # The termination guard, in numbers: every ACCEPTED shift sits between the
    # ABSOLUTE floor τ₀/N and the overlap cap τ, so the window count is bounded.
    @test all(x -> 1.0 / N - 1e-12 ≤ x ≤ 1.0 + 1e-12, diff(windowstarts(solution)))
    @test length(solution) ≤ ceil(Int, (tN - 1.0) / (1.0 / N)) + 2
    @test solution.restarts[1] == 0    # window 1 exempt
    @test any(>(0), solution.restarts) # the machinery fired
    @test all(≤(2), solution.restarts) # R caps every window
    # NOT asserted: restarts == R everywhere. A reduced shift means more
    # overlap and better seeding, so a retry can CONVERGE before exhausting R
    # — that is recovery working, not a failed guard. The old `all(==(2))`
    # assumed failure was unrecoverable, which bad_parareal's ϵ deliberately
    # is not (see its comment at the top of this file).
end

@testset "Zoom (AMoWi-3): multiplicative, per restart, on the base weight" begin
    # Unit-level direction check on the hooks:
    parareal = MoWi(bad_parareal(ψ=ψ₂, weights=Weights(w=2.0)); τ=1.0, adaptive=ZoomParameters()).parallelsolver
    state = NSDEMovingWindow.MoWiState(1.0, 1.0, 0.0, 0.0) # Zoom's hooks never touch the state
    NSDEMovingWindow.adapt_failure!(ZoomParameters(), state, parareal)
    @test parareal.tolerance.weights.w ≈ 4.0                   # ×δw⁻ = ×2: looser
    NSDEMovingWindow.adapt_success!(ZoomParameters(), state, parareal)
    @test parareal.tolerance.weights.w ≈ 1.0                   # ×δw⁺ = ×0.25: tighter
    # End-to-end: restarts fire; τ and Δτ never move; the weight has been driven.
    mowi = MoWi(bad_parareal(ψ=ψ₂, weights=Weights(w=2.0)); τ=1.0, Δτ=0.5, adaptive=ZoomParameters(R=3))
    solution = solve(problem, mowi)
    @test any(>(0), solution.restarts)
    @test all(l -> isapprox(l, 1.0; atol=1e-12), windowlengths(solution))
    @test diff(windowstarts(solution)) ≈ fill(0.5, length(solution) - 1) atol = 1e-12
    @test mowi.parallelsolver.tolerance.weights.w != 2.0
end

@testset "ensemble statistics (thesis §4.2)" begin
    # On the fixed point u ≡ 1 every window is the constant solution: the
    # ensemble mean must be exactly 1, the spread exactly 0.
    fixedpoint = Logistic(1.0, (0.0, 3.0))
    solution = solve(fixedpoint, MoWi(good_parareal(); τ=1.0, Δτ=0.5))
    @test ensemblemean(solution; samples=7) ≈ ones(7) atol = 1e-10
    @test ensemblevariance(solution; samples=7) ≈ zeros(7) atol = 1e-20
    @test ensemblesem(solution; samples=7) ≈ zeros(7) atol = 1e-20
    # General shape and the σ/√M relation on a non-trivial run:
    solution = solve(problem, MoWi(good_parareal(); τ=1.0, Δτ=0.5))
    σ² = ensemblevariance(solution, u -> u[1]; samples=11)
    @test all(≥(0), σ²)
    @test ensemblesem(solution, u -> u[1]; samples=11) ≈ sqrt.(σ² ./ length(solution))
end

@testset "evaluator edges (diagnostic use)" begin
    solution = solve(problem, MoWi(good_parareal(); τ=1.0, Δτ=0.5))
    @test solution(-5.0) == solution[1](-5.0)                  # clamped below
    @test solution(99.0) == solution[end](99.0)                # clamped above
    @test solution(0.75) == solution[2](0.75)                  # LATEST covering window wins
    @test norm(solution(1.7) - finetruth(1.7)) < 1e-6
end

@testset "container types" begin
    solution = solve(problem, MoWi(good_parareal(); τ=1.0, Δτ=0.5))
    @test eltype(solution.restarts) == Int
    for m = 1:length(solution)
        @test isconcretetype(eltype(solution[m].lastiterate.chunks)) # the hot seam stays concrete
    end
end

@testset "RecipesBase recipes (headless)" begin
    solution = solve(problem, MoWi(good_parareal(); τ=1.0, Δτ=0.5))
    for obj in (solution, NSDEBase._PhasePlot(solution), NSDEBase._Convergence(solution))
        @test RecipesBase.apply_recipe(Dict{Symbol,Any}(), obj) isa Vector{RecipesBase.RecipeData}
    end
    # The MoWi convergence recipe's 2-element-label path (K_m and R_m
    # series), exercised on purpose:
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:label => ["K", "R"]),
                                   NSDEBase._Convergence(solution)) isa Vector{RecipesBase.RecipeData}
end


end # outer testset
