"""
    solve(problem::InitialValueProblem, solver::RungeKuttaSolver) :: RungeKuttaSolution

returns the `RungeKuttaSolution` of an `InitialValueProblem`.
"""
function NSDEBase.solve(problem::InitialValueProblem, solver::RungeKuttaSolver; save_stages = false)
    solution = RungeKuttaSolution(problem, solver, save_stages)
    cache = Cache(problem, solver)
    @↓ (t₀, T) ← tspan = problem
    @↓ u, t = solution
    @↓ n = cache
    # WHILE instead of FOR loop -> adaptive methods
    N = length(t)
    while n < N && t[n] < T
        step!(solution, problem, solver, cache, save_stages)
        adaptive_step!(solution, solver, cache, save_stages)
        @↓ n = cache
    end
    # Decrease size of output vectors -> adaptive methods
    if n < N
        resize!(u, n)
        resize!(t, n)
    end
    return solution
end
