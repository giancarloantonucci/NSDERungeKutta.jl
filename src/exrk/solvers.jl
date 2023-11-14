function KassamTrefethenStabilisation(Z::AbstractMatrix, f::Function)
    N = 64
    θ = ((1:N/2) .- 0.5) * 2π/N * 1im
    α = exp.(θ)
    result = zero(Z)
    for j in eachindex(α)
        result .+= inv(α[j] * I - Z) * f(α[j]) * α[j]
    end
    return 2 * real(result)
end

@doc (@doc KassamTrefethenStabilisation) KTS(args...; kwargs...) = KassamTrefethenStabilisation(args...; kwargs...)

# Base.:*(KT::KassamTrefethenStabilisation, v::AbstractVector) = KT(v)

"""
    ETDEuler(; h::Real=0.0) :: ExplicitExponentialRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 1st-order stiff Exponential-Time-Differencing Euler method.
```
"""
function ETDEuler(Z::AbstractMatrix; h::Real=0.0)
    p = 1
    φ₁(Z) = Z \ (exp(Z) - I)
    tableau = ButcherTableau([
        0        0;
        p  KTS(φ₁);
    ])
    return ERK(tableau, h)
end

"""
    ETDRK4(; h::Real=0.0) :: ExplicitExponentialRungeKuttaSolver

returns an [`ExplicitExponentialRungeKuttaSolver`](@ref) for the 2nd-order stiff Exponential-Time-Differencing Runge-Kutta method.
"""
function ETDRK4(Z::AbstractMatrix; h::Real=0.0)
    # p = 2
    # φ₁
    # φ₂
    # φ₁
    return ExRK(h)
end
