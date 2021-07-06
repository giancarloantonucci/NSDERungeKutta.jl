"""
    AdaptiveParameters(; δ = 0.0, ϵ = 1e-5, K = 100) -> AdaptiveParameters

returns a constructor containing the parameters of an adaptive `RungeKuttaSolver`.

# Arguments
- `δ :: Real`    : absolute tolerance.
- `ϵ :: Real`    : relative tolerance.
- `K :: Integer` : maximum number of iterations.
"""
mutable struct AdaptiveParameters{δ_T, ϵ_T, K_T}
    δ::δ_T
    ϵ::ϵ_T
    K::K_T
end

AdaptiveParameters(; δ = 0.0, ϵ = 1e-5, K = 100) = AdaptiveParameters(δ, ϵ, K)

function adaptive_step!(solution, solver, cache, save_stages)
    @↓ u = solution
    @↓ h = solver
    @↓ n, m, v = cache
    if solver.adaptive isa AdaptiveParameters
        save_stages ? (@↓ k ← k[n] = solution) : (@↓ k = cache)
        @↓ b, s, p, d, q = solver.tableau
        @↓ δ, ϵ, K = solver.adaptive
        zero!(v)
        for i = 1:s
            @. v += (b[i] - d[i]) * k[i]
        end
        @. v /= δ + max(abs(u[n]), abs(u[n+1])) * ϵ
        error = norm(v) / √length(v)
        power = - 1 / (min(p, q) + 1)
        factor = 0.9 * error ^ power
        h *= max(0.5, min((m == 1 ? 2.0 : 1.0), factor))
        if error < 1 || m ≥ K
            n += 1
        else
            m += 1
        end
    else
        n += 1
    end
    @↑ cache = n, m
    @↑ solver = h
end
