function adaptivestep!(cache::ExplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, solver::ExplicitRungeKuttaSolver, adaptive::AbstractAdaptiveParameters)
    @↓ n, m, v, k = cache
    @↓ u, t = solution
    @↓ s, b, p, d, q = solver.tableau
    @↓ h, hs = solver.stepsize
    @↓ εₐ, εᵣ, Mₙ = adaptive

    save_stepsizes = !isnothing(hs)

    zero!(v)
    for i = 1:s
        @. v += (b[i] - d[i]) * k[i]
    end
    ε = εₐ + norm(u[n]) * εᵣ
    δ = hairernorm(v)
    h *= max(0.5, min(2.0, (0.35 * ε / δ) ^ (1 / (min(p, q) + 1))))

    if h ≈ zero(h)
        error("Step-size `h` too small at `t = $(t[n])`.")
    end

    if δ < ε || m ≥ Mₙ
        if save_stepsizes
            push!(hs.accepted, h)
            push!(hs.rejected, [])
        end
        m = 1
        n += 1
        if m ≥ Mₙ
            println("Maximum number of iterations reached at `t = $(t[n])`.")
        end
    else
        if save_stepsizes
            push!(hs.rejected[end], h)
        end
        m += 1
    end

    @↑ cache = n, m
    @↑ solver.stepsize = h
    return solution
end
