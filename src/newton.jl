"""
    NewtonParameters <: AbstractNewtonParameters

A composite type for an [`AbstractNewtonParameters`](@ref).

# Constructors
```julia
NewtonParameters(; ϵ=1e-3, K=10)
```

## Arguments
- `ϵ :: Real` : relative tolerance
- `K :: Integer` : maximum number of iterations.
"""
struct NewtonParameters{ϵ_T<:Real, K_T<:Integer} <: AbstractNewtonParameters
    ϵ::ϵ_T
    K::K_T
end

NewtonParameters(; ϵ=1e-3, K=10) = NewtonParameters(ϵ, K)
