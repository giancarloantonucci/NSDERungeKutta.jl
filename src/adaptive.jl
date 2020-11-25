"""
    AdaptiveParameters(δ, ϵ, K)

returns a constructor for the parameters of an adaptive method:
- `δ` : absolute tolerance.
- `ϵ` : relative tolerance.
- `K` : max number of iterations.
"""
mutable struct AdaptiveParameters{δ_T, ϵ_T, K_T}
    δ::δ_T
    ϵ::ϵ_T
    K::K_T
end
AdaptiveParameters(; δ = 0.0, ϵ = 1e-5, K = 100) = AdaptiveParameters(δ, ϵ, K)
