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
    solve!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver) :: RungeKuttaSolution

computes the `solution` of `problem` using `solver`.
"""
function NSDEBase.solve!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver)
    cache = RungeKuttaCache(problem, solver)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ u, t = solution
    @↓ n = cache
    N0 = N = length(t)
    # Integrate until reaching the end time `tN`
    while t[n] < tN
        step!(cache, solution, problem, solver)
        adaptivestep!(cache, solution, solver)
        @↓ n = cache
        # If the solution array is full and the end time hasn't been reached yet, append memory for more time steps
        if n == N && t[n] < tN
            append!(u, [similar(u[n]) for i = 1:N0])
            append!(t, similar(t, N0))
            N += N0
        end
    end
    # Resize the solution arrays to match the final number of time steps
    resize!(u, n)
    resize!(t, n)
    return solution # automatically updated
end

"""
    solve(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver) :: RungeKuttaSolution

computes the solution of `problem` using `solver`.
"""
function NSDEBase.solve(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver)
    solution = RungeKuttaSolution(problem, solver)
    solve!(solution, problem, solver)
    return solution
end
