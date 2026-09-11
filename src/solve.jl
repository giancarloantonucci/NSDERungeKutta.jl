"""
    initialize_cache(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver)

builds a cache for a Runge-Kutta solver.
"""
NSDEBase.initialize_cache(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver) = RungeKuttaCache(problem, solver)

"""
    initialize_solution(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver; kwargs...)

builds an empty solution object for a Runge-Kutta solver.
"""
NSDEBase.initialize_solution(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver; kwargs...) = RungeKuttaSolution(problem, solver; kwargs...)

function step!(cache::AbstractRungeKuttaCache, solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver)
    @↓ rhs = problem
    return step!(cache, solution, rhs, solver)
end

function adaptivestep!(cache::AbstractRungeKuttaCache, solution::AbstractRungeKuttaSolution, solver::AbstractRungeKuttaSolver)
    @↓ adaptive = solver
    return adaptivestep!(cache, solution, solver, adaptive)
end

function adaptivestep!(cache::AbstractRungeKuttaCache, solution::AbstractRungeKuttaSolution, solver::AbstractRungeKuttaSolver, adaptive::Nothing)
    cache.n += 1
    return solution
end

"""
    solve!(cache::AbstractRungeKuttaCache, solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver) :: RungeKuttaSolution

computes the `solution` of `problem` using `solver` and a pre-allocated `cache`.
"""
function NSDEBase.solve!(cache::AbstractRungeKuttaCache, solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver)
    # Reset cache state for reuse
    cache.n = 1
    cache.m = 1
    cache.e[] = 0.0

    @↓ u0, (t0, tN) ← tspan = problem
    @↓ u, t, k = solution
    N0 = N = length(t)

    # Reused solutions must not keep the previous run's stage history:
    if !(k isa Nothing)
        empty!(k)
    end

    # CRITICAL FIX: Inject the new chunk's initial conditions into the reused arrays
    copyto!(u[1], u0)
    t[1] = t0

    # Number of steps after which the fixed-step schedule ends on `tN`, or 0
    # when the step does not divide the span (or the solver is adaptive, whose
    # schedule is not known in advance):
    M_last = fixedstep_count(solver, t0, tN)

    # Integrate until reaching the end time `tN`
    while t[cache.n] < tN
        step!(cache, solution, problem, solver)

        n_before = cache.n
        adaptivestep!(cache, solution, solver)

        # Fixed-step schedule: the step that the schedule says is the last one
        # lands on `tN` exactly. The time sum can come to rest an ulp either
        # side of `tN`; below it the loop would take a whole extra step, above
        # it a lookup at `tN` would interpolate rather than return the last
        # node. The decision rests on the STEP COUNT, never on a time
        # tolerance alone: a few ulps of absolute time can exceed the step
        # itself far from the origin, and would then snap a step that is not
        # the last one.
        if M_last > 0 && cache.n == M_last + 1 && t[cache.n] != tN
            t[cache.n] = tN
        end

        # Save stages for dense output, but only for ACCEPTED steps: pushing on
        # every attempt lets rejected stages drift the history out of line with
        # the accepted grid.
        if !(k isa Nothing) && cache.n > n_before
            push!(k, copy.(cache.k))
        end

        if cache.n == N && t[cache.n] < tN
            append!(u, [similar(u[cache.n]) for i = 1:N0])
            append!(t, similar(t, N0))

            if !(k isa Nothing)
                sizehint!(k, length(k) + N0)
            end
            N += N0
        end
    end

    resize!(u, cache.n)
    resize!(t, cache.n)
    if !(k isa Nothing)
        resize!(k, cache.n - 1) # stages live on intervals: one fewer than grid points
    end
    return solution
end

"""
    fixedstep_count(solver, t0, tN) :: Int

the number of steps `M` after which a FIXED-step schedule `t0 + M h` reaches
`tN`, or `0` when it does not: the step does not divide the span, or the
solver is adaptive (its schedule is not known in advance).

Divisibility is judged on the LOCAL quotient `q = (tN − t0) / h` alone: `q`
must sit within a few ulps of itself of an integer. Nothing here depends on
the absolute time origin — an ulp of `|t0|` can be larger than `h` far from
the origin, and any allowance in those units would wave through a schedule
that misses `tN` by a real fraction of a step. Where the span subtraction is
inexact the test can only fail to recognise a dividing schedule, which falls
back to the ordinary end-of-loop policy (a final node at or beyond `tN`).
"""
function fixedstep_count(solver::AbstractRungeKuttaSolver, t0::Real, tN::Real)
    hasproperty(solver, :adaptive) && !(solver.adaptive isa Nothing) && return 0
    h = solver.stepsize.h
    h > 0 || return 0
    q = (float(tN) - float(t0)) / float(h)
    isfinite(q) && 1 ≤ q ≤ typemax(Int) ÷ 2 || return 0
    M = round(Int, q)
    abs(q - M) ≤ 8 * eps(max(q, one(q))) || return 0
    return M
end

"""
    solve!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver) :: RungeKuttaSolution

computes the `solution` of `problem` using `solver`, allocating a new cache.
"""
function NSDEBase.solve!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver)
    cache = RungeKuttaCache(problem, solver)
    return NSDEBase.solve!(cache, solution, problem, solver)
end

"""
    solve(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver; dense::Bool=false, kwargs...) :: RungeKuttaSolution

computes the solution of `problem` using `solver`.
"""
function NSDEBase.solve(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver; dense::Bool=false)
    solution = RungeKuttaSolution(problem, solver; dense=dense)
    NSDEBase.solve!(solution, problem, solver)
    return solution
end
