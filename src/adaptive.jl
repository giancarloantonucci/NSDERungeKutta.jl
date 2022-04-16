"""
    AdaptiveParameters <: AbstractAdaptiveParameters

A composite type for the parameters of an adaptive [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
AdaptiveParameters(; atol::Real=0.0, rtol::Real=1e-5, nits::Integer=100)
```

## Arguments
- `atol :: Real` : absolute tolerance.
- `rtol :: Real` : relative tolerance.
- `nits :: Integer` : maximum number of iterations.
"""
struct AdaptiveParameters{atol_T<:Real, rtol_T<:Real, nits_T<:Integer} <: AbstractAdaptiveParameters
    atol::atol_T
    rtol::rtol_T
    nits::nits_T
end
AdaptiveParameters(; atol::Real=0.0, rtol::Real=1e-5, nits::Integer=100) = AdaptiveParameters(atol, rtol, nits)
