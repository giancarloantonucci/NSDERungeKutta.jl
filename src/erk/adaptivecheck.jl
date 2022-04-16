function adaptivecheck!(cache::ExplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, solver::ExplicitRungeKuttaSolver, adaptive::AbstractAdaptiveParameters)
    @↓ n, m, v, k = cache
    @↓ u = solution
    @↓ tableau, stepsize = solver
    @↓ s, b, p, d, q = tableau
    @↓ h = stepsize
    @↓ atol, rtol, nits = adaptive
    zero!(v)
    for i = 1:s
        @. v += (b[i] - d[i]) * k[i]
    end
    tol = atol + norm(u[n]) * rtol
    err = hairernorm(v)
    pwr = 1 / (min(p, q) + 1)
    r = max(0.5, min(2.0, (0.35 * tol / err) ^ pwr))
    h *= r
    if err < tol || m ≥ nits
        n += 1
    else
        m += 1
    end
    @↑ cache = n, m
    @↑ stepsize = h
    return solution
end
