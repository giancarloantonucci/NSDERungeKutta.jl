"""
    initialize_cache(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver)

builds a cache for a Runge-Kutta solver.
"""
function NSDEBase.initialize_cache(
    problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver
)
    return RungeKuttaCache(problem, solver)
end

"""
    initialize_solution(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver; kwargs...)

builds an empty solution object for a Runge-Kutta solver.
"""
function NSDEBase.initialize_solution(
    problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver; kwargs...
)
    return RungeKuttaSolution(problem, solver; kwargs...)
end

@inline function step!(
    cache::AbstractRungeKuttaCache,
    solution::AbstractRungeKuttaSolution,
    problem::AbstractInitialValueProblem,
    solver::AbstractRungeKuttaSolver,
    h::Real=solver.stepsize.h,
)
    @↓ rhs = problem
    return step!(cache, solution, rhs, solver, h)
end

function adaptivestep!(
    cache::AbstractRungeKuttaCache,
    solution::AbstractRungeKuttaSolution,
    solver::AbstractRungeKuttaSolver,
)
    @↓ adaptive = solver
    return adaptivestep!(cache, solution, solver, adaptive)
end

function adaptivestep!(
    cache::AbstractRungeKuttaCache,
    solution::AbstractRungeKuttaSolution,
    solver::AbstractRungeKuttaSolver,
    adaptive::Nothing,
)
    cache.n += 1
    return solution
end

"""
    solve!(cache::AbstractRungeKuttaCache, solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver) :: RungeKuttaSolution

computes the `solution` of `problem` using `solver` and a pre-allocated `cache`.
"""
function NSDEBase.solve!(
    cache::AbstractRungeKuttaCache,
    solution::AbstractRungeKuttaSolution,
    problem::AbstractInitialValueProblem,
    solver::AbstractRungeKuttaSolver,
)
    check_dense_support(solver, solution.k !== nothing)
    if solver.adaptive === nothing
        return solve_fixed!(
            cache,
            solution,
            problem,
            solver,
            StepSchedule(problem.tspan[2] - problem.tspan[1], solver.stepsize.h),
        )
    end
    reset_solution!(cache, solution, problem)
    tN = problem.tspan[2]
    h = solver.stepsize.h
    while solution.t[cache.n] < tN
        n = cache.n
        ensure_capacity!(solution, n + 1, problem.u0)
        remaining = tN - solution.t[n]
        rounding = max(8eps(float(tN)), 1e-6 * h)
        finishing = remaining <= h + rounding
        hused = finishing ? remaining : h
        solution.t[n] + typeof(solution.t[n])(hused) > solution.t[n] ||
            throw(StepFailure("adaptive step cannot advance time at this precision"))
        compensation = cache.e[]
        step!(cache, solution, problem, solver, hused)
        h = adaptivestep!(cache, solution, solver, solver.adaptive, hused)
        if cache.n > n
            finishing && (solution.t[cache.n] = tN)
            solution.k === nothing || push!(solution.k, copy.(cache.k))
        else
            cache.e[] = compensation
        end
    end
    return finish_solution!(solution, cache.n)
end

"""
    solve!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver) :: RungeKuttaSolution

computes the `solution` of `problem` using `solver`, allocating a new cache.
"""
function NSDEBase.solve!(
    solution::AbstractRungeKuttaSolution,
    problem::AbstractInitialValueProblem,
    solver::AbstractRungeKuttaSolver,
)
    cache = RungeKuttaCache(problem, solver)
    return NSDEBase.solve!(cache, solution, problem, solver)
end

"""
    solve(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver; dense::Bool=false, kwargs...) :: RungeKuttaSolution

computes the solution of `problem` using `solver`. `dense=true` stores stage
history for dense interpolation. IMEX solvers currently require `dense=false`;
their returned solutions still support linear interpolation.
"""
function NSDEBase.solve(
    problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver; dense::Bool=false
)
    solution = RungeKuttaSolution(problem, solver; dense=dense)
    NSDEBase.solve!(solution, problem, solver)
    return solution
end
