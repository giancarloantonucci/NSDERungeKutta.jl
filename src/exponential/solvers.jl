"""
    ExponentialRK4(; h::Real=0.0) :: ExplicitExponentialRungeKuttaSolver
    ERK4(args...; kwargs...)

returns an [`ExplicitExponentialRungeKuttaSolver`](@ref) for the 4th-order EERK method.
"""
ExponentialRK4(; h::Real=0.0) = EERK(h)
@doc (@doc ExponentialRK4) ERK4(args...; kwargs...) = ExponentialRK4(args...; kwargs...)
