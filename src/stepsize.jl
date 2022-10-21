"""
    StepSize <: AbstractStepSize

A composite type for the step-size a Runge-Kutta solver.

# Constructors
```julia
StepSize(h::Real)
```

# Functions
[`stepsize`](@ref) : returns step-size.
"""
mutable struct StepSize{h_T<:Real} <: AbstractStepSize
    h::h_T
end

#---------------------------------- FUNCTIONS ----------------------------------

"""
    stepsize(solver::AbstractRungeKuttaSolver) :: Real

returns the step-size of a `solver`.
"""
stepsize(solver::AbstractRungeKuttaSolver) = solver.stepsize.h
