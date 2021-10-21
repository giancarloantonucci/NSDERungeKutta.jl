"""
    ExplicitExponentialRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for an exponential explicit [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
ExplicitExponentialRungeKuttaSolver(h)
EERK(args...; kwargs...)
```

# Arguments
- `h :: Real` : step-size.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
mutable struct ExplicitExponentialRungeKuttaSolver{stepsize_T} <: AbstractRungeKuttaSolver
    stepsize::stepsize_T
end

function ExplicitExponentialRungeKuttaSolver(h::Real)
    return ExplicitExponentialRungeKuttaSolver(StepSize(h))
end

@doc (@doc ExplicitExponentialRungeKuttaSolver) function EERK(args...; kwargs...)
    return ExplicitExponentialRungeKuttaSolver(args...; kwargs...)
end

#####
##### Methods
#####

# POSSIBLE SPECIALISATION
# function (solver::ExplicitExponentialRungeKuttaSolver)(problem::AbstractInitialValueProblem)
#     ...
# end
