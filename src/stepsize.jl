# TODO: StepSizes → StartStepSizes

mutable struct StepSizes{accepted_T<:AbstractVector{<:Real}, rejected_T<:AbstractVector{<:AbstractVector{<:Real}}}
    accepted::accepted_T
    rejected::rejected_T
end

function StepSizes(h::Real)
    accepted = typeof(h)[]
    push!(accepted, h)
    rejected = typeof(accepted)[]
    push!(rejected, [])
    return StepSizes(accepted, rejected)
end

"""
    StepSize <: AbstractStepSize

A composite type for the step-size a Runge-Kutta solver.

# Constructors
```julia
StepSize(h::Real)
```

# Functions
[`stepsize`](@ref) : returns (last) step-size
"""
mutable struct StepSize{h_T<:Real, hs_T<:Union{StepSizes,Nothing}} <: AbstractStepSize
    h::h_T
    hs::hs_T
end

function StepSize(h::Real; save_stepsizes::Bool=false)
    if save_stepsizes
        return StepSize(h, StepSizes(h))
    else
        return StepSize(h, nothing)
    end
end

#---------------------------------- FUNCTIONS ----------------------------------

"""
    stepsize(solver::AbstractRungeKuttaSolver) :: Real

returns the step-size of a `solver`.
"""
stepsize(solver::AbstractRungeKuttaSolver) = solver.stepsize.h
