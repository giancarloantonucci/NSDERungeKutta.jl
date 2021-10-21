"""
    ExponentialRK4(; h::Real=0.0) :: ExplicitExponentialRungeKuttaSolver

returns an [`ExplicitExponentialRungeKuttaSolver`](@ref) for the 4th-order EERK method.
"""
function ExponentialRK4(; h::Real=0.0)
    return EERK(h)
end

@doc (@doc ExponentialRK4) function ERK4(args...; kwargs...)
    return ExponentialRK4(args...; kwargs...)
end
