"""
    ExplicitExponentialRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for an exponential explicit [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
ExplicitExponentialRungeKuttaSolver(stepsize::StepSize)
EERK(args...; kwargs...)
```

# Methods

    (solver::ExplicitExponentialRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    
returns the `solution` of a `problem` using `solver`.
    
    (solver::ExplicitExponentialRungeKuttaSolver)(problem::AbstractInitialValueProblem) :: RungeKuttaSolution

returns the [`RungeKuttaSolution`](@ref) of a `problem` using `solver`.
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

(solver::ExplicitExponentialRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) = solve!(solution, problem, solver)
(solver::ExplicitExponentialRungeKuttaSolver)(problem::AbstractInitialValueProblem) = solve(problem, solver)
