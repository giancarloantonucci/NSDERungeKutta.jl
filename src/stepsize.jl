"""
    StepSize <: AbstractStepSize

A composite type for an [`AbstractStepSize`](@ref).

# Constructors
```julia
StepSize(; h::Real)
```

# Functions
[`stepsize`](@ref) : returns step-size.
"""
mutable struct StepSize{h_T<:Real} <: AbstractStepSize
    h::h_T
end

"""
    stepsize(solver::AbstractRungeKuttaSolver)

returns the step-size of a `solver`.
"""
stepsize(solver::AbstractRungeKuttaSolver) = solver.stepsize.h
