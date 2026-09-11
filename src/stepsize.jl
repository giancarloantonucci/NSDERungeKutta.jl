# NSDERungeKutta/src/stepsize.jl

# TODO: StepSizes → StartStepSizes

"""
    StepSizes

The REALISED step-size history of an adaptive run (opt-in via
`save_stepsizes = true`):

- `accepted[i]` is exactly the `i`-th step taken, so
  `hs.accepted == diff(solution.t)` — the record is what happened, never the
  controller's proposal for the next attempt.
- `rejected[i]` holds the step sizes tried and REJECTED before accepted step
  `i` (the values that failed, not the reduced retries they triggered).
  A trailing empty bin is opened at every acceptance, so
  `length(rejected) == length(accepted) + 1` and the last bin is empty at the
  end of a completed run.
"""
mutable struct StepSizes{accepted_T<:AbstractVector{<:Real}, rejected_T<:AbstractVector{<:AbstractVector{<:Real}}}
    accepted::accepted_T
    rejected::rejected_T
end

function StepSizes(h::Real)
    # `accepted` starts EMPTY: pre-seeding it with the configured h recorded a
    # step that was never necessarily taken (the first attempt can be
    # rejected). `rejected` starts with one empty bin — the container for
    # failures before the first acceptance.
    accepted = typeof(h)[]
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
