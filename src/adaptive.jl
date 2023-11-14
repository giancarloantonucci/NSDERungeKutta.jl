"""
    AdaptiveParameters <: AbstractAdaptiveParameters

A composite type for the parameters of an adaptive [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
AdaptiveParameters(εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100)
```

## Arguments
- `εₐ :: Real` : absolute tolerance
- `εᵣ :: Real` : relative tolerance
- `Mₙ :: Integer` : maximum number of iterations
"""
struct AdaptiveParameters{εₐ_T<:Real, εᵣ_T<:Real, Mₙ_T<:Integer} <: AbstractAdaptiveParameters
    εₐ :: εₐ_T
    εᵣ :: εᵣ_T
    Mₙ :: Mₙ_T
end

AdaptiveParameters(; εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100) = AdaptiveParameters(εₐ, εᵣ, Mₙ)
