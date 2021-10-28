"""
    StepSize <: AbstractStepSize

A composite type for an [`AbstractStepSize`](@ref).

# Constructors
```julia
StepSize(; h)
```

# Arguments
- `h :: Real` : step-size.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
mutable struct StepSize{h_T} <: AbstractStepSize
    h::h_T
end

stepsize(solver::AbstractRungeKuttaSolver) = solver.stepsize.h