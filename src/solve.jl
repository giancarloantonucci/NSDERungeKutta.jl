function check_adaptive!(cache, solution, solver)
    if solver.adaptive isa AdaptiveParameters
        @↓ n, m, k = cache
        @↓ u = solution
        @↓ tableau, p, p̂, h, adaptive = solver
        @↓ b, b̂, s = tableau
        @↓ δ, ϵ, K = adaptive
        v = copy(u[Block(n+1)]) # to avoid allocs
        L = length(v)
        zero!(v)
        for i = 1:s
            @. v += (b[i] - b̂[i]) * k[Block(i)]
        end
        @. v /= δ + abs(u[Block(n)]) * ϵ
        err = sqrt(sum(v.^2) / L)
        q = 1 / (min(p, p̂) + 1)
        h *= max(0.5, min((m == 1 ? 2.0 : 1.0), 0.9 / (err ^ q)))
        @↑ solver = h
        if err < 1.0 || m ≥ K
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
    @↓ tspan = problem
    tN = tspan[2]
    solution = RungeKuttaSolution(problem, solver)
    @↓ u, t = solution
    N = length(t)
    cache = RungeKuttaCache(problem, solver)
    @↓ n = cache
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
