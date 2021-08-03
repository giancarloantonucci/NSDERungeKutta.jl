"""
    StepSize

A composite type for the step-size of a [`RungeKuttaSolver`](@ref).

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
mutable struct StepSize{h_T}
    h::h_T
end

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

"""
    show(io::IO, stepsize::StepSize)

prints a full description of `stepsize` and its contents to a stream `io`.
"""
Base.show(io::IO, stepsize::StepSize) = _show(io, stepsize)

"""
    summary(io::IO, stepsize::StepSize)

prints a brief description of `stepsize` to a stream `io`.
"""
Base.summary(io::IO, stepsize::StepSize) = _summary(io, stepsize)
