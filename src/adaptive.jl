"""
    AdaptiveParameters(; δ = 0.0, ϵ = 1e-5, K = 100) -> AdaptiveParameters

returns a constructor for the parameters of an adaptive `RungeKuttaSolver`.

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
