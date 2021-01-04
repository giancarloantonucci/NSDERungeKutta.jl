function check_adaptive!(cache, solution, solver)
    if solver.adaptive isa AdaptiveParameters
        @↓ n, m, v, k = cache
        @↓ u = solution
        @↓ tableau, adaptive = solver
        @↓ b, s, p, d, q = tableau
        @↓ δ, ϵ, K = adaptive
        L = length(v)
        zero!(v)
        for i = 1:s
            @. v += (b[i] - d[i]) * k[Block(i)]
        end
        @. v /= δ + max(abs(u[Block(n)]), abs(u[Block(n+1)])) * ϵ
        error = norm(v) / sqrt(L)
        power = -1 / (min(p, q) + 1)
        solver.h *= max(0.5, min((m == 1 ? 2.0 : 1.0), 0.9 * error ^ power))
        if error < 1 || m ≥ K
            cache.n += 1
        else
            cache.m += 1
        end
    else
        cache.n += 1
    end
end

"""
    solve(problem::InitialValueProblem, solver::RungeKuttaSolver) -> RungeKuttaSolution

returns the `RungeKuttaSolution` of an `InitialValueProblem`.
"""
function NSDEBase.solve(problem::InitialValueProblem, solver::RungeKuttaSolver)
    solution = RungeKuttaSolution(problem, solver)
    cache = RungeKuttaCache(problem, solver)
    @↓ (t0, tN) ← tspan = problem
    @↓ u, t = solution
    @↓ n = cache
    N = length(t)
    while n < N && t[n] < tN
        step!(cache, solution, problem, solver)
        check_adaptive!(cache, solution, solver)
        @↓ n = cache
    end
    if n < N
        resize!(u.blocks, n)
        resize!(t, n)
    end
    return solution
end
