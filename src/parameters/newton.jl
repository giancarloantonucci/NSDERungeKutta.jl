"""
    NewtonParameters(; ϵ = 1e-3, K = 10) <: RungeKuttaParameters

returns a constructor containing the parameters of simplified Newton in `ImplicitRungeKuttaSolver`.

# Arguments
- `ϵ :: Real`    : relative tolerance
- `K :: Integer` : maximum number of iterations.
"""
mutable struct NewtonParameters{ϵ_T, K_T} <: RungeKuttaParameters
    ϵ::ϵ_T
    K::K_T
end

NewtonParameters(; ϵ = 1e-3, K = 10) = NewtonParameters(ϵ, K)
