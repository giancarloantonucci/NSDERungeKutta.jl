function step!(cache::AbstractRungeKuttaCache, solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver)
    @↓ rhs = problem
    return step!(cache, solution, rhs, solver)
end

function adaptivecheck!(cache::AbstractRungeKuttaCache, solution::AbstractRungeKuttaSolution, solver::AbstractRungeKuttaSolver)
    @↓ adaptive = solver
    return adaptivecheck!(cache, solution, solver, adaptive)
end

function adaptivecheck!(cache::AbstractRungeKuttaCache, solution::AbstractRungeKuttaSolution, solver::AbstractRungeKuttaSolver, adaptive::Nothing)
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
    N = N₀ = length(t)
    while t[n] < tN
        step!(cache, solution, problem, solver)
        adaptivecheck!(cache, solution, solver)
        @↓ n = cache
        @↓ h = solver.stepsize
        if t[n] ≈ tN || h ≈ zero(h)
            break
        end
        if n == N && t[n] < tN
            append!(u, similar(u, N₀))
            append!(t, similar(t, N₀))
            N += N₀
        end
    end
    if n < N
        resize!(u, n)
        resize!(t, n)
    end
    return solution
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
