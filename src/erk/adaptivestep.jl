function adaptivestep!(cache::ExplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, solver::ExplicitRungeKuttaSolver, adaptive::AbstractAdaptiveParameters)
    @↓ n, m, v, k = cache
    @↓ u, t = solution
    @↓ tableau, stepsize = solver
    @↓ s, b, p, d, q = tableau
    @↓ h = stepsize
    @↓ ϵₐ ← atol, ϵᵣ ← rtol, M ← nits = adaptive
    zero!(v)
    for i = 1:s
        @. v += (b[i] - d[i]) * k[i]
    end
    ϵ = ϵₐ + norm(u[n]) * ϵᵣ
    err = hairernorm(v)
    h *= max(0.5, min(2.0, (0.35 * ϵ / err) ^ (1 / (min(p, q) + 1))))
    if h ≈ zero(h)
        error("Stepsize `h` has become too small at `t = $(t[n])`.")
    end
    if err < ϵ || m ≥ M
        m = 1
        n += 1
        if m ≥ M
            println("Adaptive loop has reached `nits` at `t = $(t[n])`.")
        end
    else
        m += 1
    end
    @↑ cache = n, m
    @↑ stepsize = h
    return solution
end
