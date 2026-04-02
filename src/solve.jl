# NSDERungeKutta/src/solve.jl

# FILE: /Users/antonucci/.julia/dev/NSDERungeKutta/src/solve.jl

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

    # CRITICAL FIX: Inject the new chunk's initial conditions into the reused arrays
    copyto!(u[1], u0)
    t[1] = t0

    # Integrate until reaching the end time `tN`
    while t[cache.n] < tN
        step!(cache, solution, problem, solver)

        # Save stages for dense output
        if !(k isa Nothing)
            push!(k, copy.(cache.k))
        end

        adaptivestep!(cache, solution, solver)

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
    return solution
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
function NSDEBase.solve(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver; dense::Bool=false, kwargs...)
    solution = RungeKuttaSolution(problem, solver; dense=dense)
    NSDEBase.solve!(solution, problem, solver; kwargs...)
    return solution
end
