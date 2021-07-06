"""
    solve(problem::InitialValueProblem, solver::RungeKuttaSolver) -> RungeKuttaSolution

returns the `RungeKuttaSolution` of an `InitialValueProblem`.
"""
function NSDEBase.solve(problem::InitialValueProblem, solver::RungeKuttaSolver; save_stages = false)
    solution = RungeKuttaSolution(problem, solver, save_stages)
    cache = Cache(problem, solver)
    @↓ (t0, tN) ← tspan = problem
    @↓ u, t = solution
    @↓ n = cache
    N = length(t)
    while n < N && t[n] < tN
        step!(solution, problem, solver, cache, save_stages)
        adaptive_step!(solution, solver, cache, save_stages)
        @↓ n = cache
    end
    if n < N
        resize!(u, n)
        resize!(t, n)
    end
    return solution
end
