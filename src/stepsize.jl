"A run's accepted step sizes and rejected attempts, stored on the solution."
struct StepSizes{T<:Real}
    accepted::Vector{T}
    rejected::Vector{Vector{T}}
end
StepSizes(h::Real) = StepSizes(typeof(h)[], [typeof(h)[]])

"Immutable initial step-size settings. Realised histories belong to the solution."
struct StepSize{T<:Real} <: AbstractStepSize
    h::T
    save_stepsizes::Bool
    function StepSize(h::Real; save_stepsizes::Bool=false)
        isfinite(h) && h > 0 || throw(ArgumentError("step size must be finite and positive"))
        return new{typeof(float(h))}(float(h), save_stepsizes)
    end
end

"""
    stepsize(solver::AbstractRungeKuttaSolver) :: Real

returns the configured (initial) step size of `solver`.
"""
stepsize(solver::AbstractRungeKuttaSolver) = solver.stepsize.h
