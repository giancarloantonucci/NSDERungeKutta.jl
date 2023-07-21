function adaptivestep!(cache::ExplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, solver::ExplicitRungeKuttaSolver, adaptive::AbstractAdaptiveParameters)
    @↓ n, m, v, k = cache
    @↓ u, t = solution
    @↓ tableau, stepsize = solver
    @↓ s, b, p, d, q = tableau
    @↓ h = stepsize
    @↓ atol, rtol, nits = adaptive
    # Compute the error estimate and target error
    zero!(v)
    for i = 1:s
        @. v += (b[i] - d[i]) * k[i]
    end
    ϵ = atol + norm(u[n]) * rtol
    err = hairernorm(v)
    # Adjust the step-size based on the error and the order of the method
    h *= max(0.5, min(2.0, (0.35 * ϵ / err) ^ (1 / (min(p, q) + 1))))
    if h ≈ zero(h)
        error("Step-size `h` too small at `t = $(t[n])`.")
    end
    # Check if the error is within the tolerance or if the maximum iterations reached
    if err < ϵ || m ≥ nits
        m = 1
        n += 1
        if m ≥ nits
            println("Maximum number of iterations `nits` hit at `t = $(t[n])`.")
        end
    else
        m += 1
    end
    # Update the stepsize and the cache
    @↑ cache = n, m
    @↑ stepsize = h
    return solution
end
