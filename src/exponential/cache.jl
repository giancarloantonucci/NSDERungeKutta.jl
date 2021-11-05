"""
    ExplicitExponentialRungeKuttaCache <: AbstractRungeKuttaCache

A composite type for the [`AbstractRungeKuttaCache`](@ref) of an [`ExplicitExponentialRungeKuttaSolver`](@ref).

# Constructors
```julia
ExplicitExponentialRungeKuttaCache(n, Q, α, β, γ, E, E2)
ExplicitExponentialRungeKuttaCache(problem, solver)
```

## Arguments
- `problem :: AbstractInitialValueProblem`
- `solver :: ExplicitRungeKuttaSolver`

# Functions
- [`RungeKuttaCache`](@ref) : alternative caller.
"""
mutable struct ExplicitExponentialRungeKuttaCache{n_T, Q_T, α_T, β_T, γ_T, E_T, E2_T} <: AbstractRungeKuttaCache
    n::n_T
    Q::Q_T
    α::α_T
    β::β_T
    γ::γ_T
    E::E_T
    E2::E2_T
end

function ExplicitExponentialRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver)
    @↓ L = problem.rhs
    @↓ h = solver.stepsize
    #
    n = 1
    #
    M = 16
    r = @. exp(π * 1im * ((1:M) - 0.5) / M)
    #
    Z = h * L .+ r'
    Q = vec(h / M * real.(sum(@. (exp(Z / 2) - 1.0) / Z; dims=2)))
    α = vec(h / M * real.(sum(@. (-4.0 - Z + exp(Z) * (4.0 - 3Z + Z^2)) / Z^3; dims=2)))
    β = vec(h / M * real.(sum(@. (2.0 + Z + exp(Z) * (-2.0 + Z)) / Z^3; dims=2)))
    γ = vec(h / M * real.(sum(@. (-4.0 - 3Z - Z^2 + exp(Z) * (4.0 - Z)) / Z^3; dims=2)))
    # L ≡ diagonal matrix
    E = @. exp(h * L)
    E2 = @. exp(h / 2 * L)
    return ExplicitExponentialRungeKuttaCache(n, Q, α, β, γ, E, E2)
end

#####
##### Functions
#####

RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver) = ExplicitExponentialRungeKuttaCache(problem, solver)
