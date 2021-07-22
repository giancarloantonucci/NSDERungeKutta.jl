"""
    StepSize(; h) <: RungeKuttaParameters

returns a constructor containing the step-size of a `RungeKuttaSolver`.

# Arguments
- `h :: Real` : step-size.
"""
mutable struct StepSize{h_T} <: RungeKuttaParameters
    h::h_T
end
