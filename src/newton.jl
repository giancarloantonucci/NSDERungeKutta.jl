"""
    NewtonParameters <: AbstractNewtonParameters

A composite type for the parameters of simplified Newton.

# Constructors
```julia
NewtonParameters(; εᵣ=1e-3, Mₙ=10)
```

## Arguments
- `εᵣ :: Real` : relative tolerance
- `Mₙ :: Integer` : maximum number of iterations
"""
struct NewtonParameters{εᵣ_T<:Real, Mₙ_T<:Integer} <: AbstractNewtonParameters
    εᵣ::εᵣ_T
    Mₙ::Mₙ_T
end

NewtonParameters(; εᵣ::Real=1e-3, Mₙ::Integer=10) = NewtonParameters(εᵣ, Mₙ)
