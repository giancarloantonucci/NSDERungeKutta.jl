"""
    NewtonParameters <: AbstractNewtonParameters

A composite type for an [`AbstractNewtonParameters`](@ref).

# Constructors
```julia
NewtonParameters(; ϵ = 1e-3, K = 10)
```

# Arguments
- `ϵ :: Real` : relative tolerance
- `K :: Integer` : maximum number of iterations.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct NewtonParameters{ϵ_T, K_T} <: AbstractNewtonParameters
    ϵ::ϵ_T
    K::K_T
end

NewtonParameters(; ϵ=1e-3, K=10) = NewtonParameters(ϵ, K)
