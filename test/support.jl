# Shared fixtures and helpers for this package test suite.
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
    hs = [h0 / 2.0^j for j in 0:(refinements - 1)]
    errs = map(hs) do h
        solution = solve(problem, make_solver(h))
        return norm(solution.u[end] - exact(solution.t[end])) + 1e-300
    end
    x, y = log.(hs), log.(errs)
    x̄, ȳ = sum(x) / length(x), sum(y) / length(y)
    lsq = sum((x .- x̄) .* (y .- ȳ)) / sum(abs2, x .- x̄)
    j = length(errs)
    if errs[j] < 50 * eps() && length(errs) > 2
        j -= 1 # finest error is machine noise: use the previous rung
    end
    ratio = log(errs[j - 1] / errs[j]) / log(hs[j - 1] / hs[j])
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

erk_fixed = [
    Euler,
    Heun2,
    Midpoint,
    Ralston2,
    Heun3,
    RungeKutta3,
    Ralston3,
    SSPRK3,
    Ralston4,
    RungeKutta4,
    Rule38,
    Butcher5,
    KuttaNystrom5,
    Butcher6,
    Butcher7,
]
erk_embedded = [HeunEuler, BogackiShampine, Fehlberg45, DormandPrince54, Verner65, Fehlberg78]
dirk_zoo = [
    BackwardEuler,
    ImplicitMidpoint,
    SDIRK2,
    LobattoIII2,
    CrankNicolson,
    SDIRK3,
    RadauI3,
    RadauII3,
    SDIRK4,
    LobattoIII4,
]
irk_zoo = [
    LobattoIIIC2,
    RadauIA3,
    RadauIIA3,
    GaussLegendre4,
    LobattoIIIA4,
    LobattoIIIB4,
    LobattoIIIC4,
    RadauI5,
    RadauIA5,
    RadauII5,
    RadauIIA5,
    GaussLegendre6,
]
ierk_zoo = [IMEXEuler, IMEXSSP2_222, IMEXSSP2_322, IMEXSSP2_332, IMEXSSP3_332]
