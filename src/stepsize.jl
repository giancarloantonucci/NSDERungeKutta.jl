"""
    StepSize(; h) -> StepSize

returns a constructor containing the step-size of a `RungeKuttaSolver`.

# Arguments
- `h :: Real` : step-size.
"""
mutable struct StepSize{h_T}
    h::h_T
end
