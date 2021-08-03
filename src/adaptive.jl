"""
    AdaptiveParameters

A composite type for the parameters of an adaptive [`RungeKuttaSolver`](@ref).

# Constructors
```julia
AdaptiveParameters(; δ = 0.0, ϵ = 1e-5, K = 100)
```

# Arguments
- `δ :: Real` : absolute tolerance.
- `ϵ :: Real` : relative tolerance.
- `K :: Integer` : maximum number of iterations.

# Methods
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
mutable struct AdaptiveParameters{δ_T, ϵ_T, K_T}
    δ::δ_T
    ϵ::ϵ_T
    K::K_T
end

AdaptiveParameters(; δ = 0.0, ϵ = 1e-5, K = 100) = AdaptiveParameters(δ, ϵ, K)

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

function adaptive_step!(solution::RungeKuttaSolution, solver::RungeKuttaSolver, adaptive::AdaptiveParameters, cache::RungeKuttaCache)
    @↓ u = solution
    @↓ h = solver.stepsize
    @↓ n, m, v = cache
    k = solution.k isa Nothing ? cache.k : solution.k[n]
    @↓ s, b, p, d, q = solver.tableau
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
    @↑ cache = n, m
    @↑ solver.stepsize = h
    return solution
end

function adaptive_step!(solution::RungeKuttaSolution, solver::RungeKuttaSolver, adaptive::Nothing, cache::RungeKuttaCache)
    cache.n += 1
    return solution
end

adaptive_step!(solution::RungeKuttaSolution, solver::RungeKuttaSolver, cache::RungeKuttaCache) = adaptive_step!(solution, solver, solver.adaptive, cache)

"""
    show(io::IO, adaptive::AdaptiveParameters)

prints a full description of `adaptive` and its contents to a stream `io`.
"""
Base.show(io::IO, adaptive::AdaptiveParameters) = _show(io, adaptive)

"""
    summary(io::IO, adaptive::AdaptiveParameters)

prints a brief description of `adaptive` to a stream `io`.
"""
Base.summary(io::IO, adaptive::AdaptiveParameters) = _summary(io, adaptive)
