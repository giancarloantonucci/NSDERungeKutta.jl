"""
    ExplicitExponentialRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for an exponential explicit [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
ExplicitExponentialRungeKuttaSolver(tableau, stepsize[, adaptive])
EERK(args...; kwargs...)
```

# Arguments
- `tableau :: AbstractButcherTableau`
- `stepsize :: Union{AbstractStepSize, Real}`
- `adaptive :: Union{AbstractAdaptiveParameters, Nothing}`

# Methods

    (solver::ExplicitExponentialRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    (solver::ExplicitExponentialRungeKuttaSolver)(problem::AbstractInitialValueProblem) :: RungeKuttaSolution

returns the `solution` of a `problem` using `solver`.
"""
struct ExplicitExponentialRungeKuttaSolver{tableau_T<:AbstractButcherTableau, stepsize_T<:AbstractStepSize, adaptive_T<:Union{AbstractAdaptiveParameters, Nothing}} <: AbstractRungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    adaptive::adaptive_T
end

ExplicitExponentialRungeKuttaSolver(tableau::AbstractButcherTableau, h::Real, adaptive::Union{AbstractAdaptiveParameters, Nothing}) = ExplicitExponentialRungeKuttaSolver(tableau, StepSize(h), adaptive)
ExplicitExponentialRungeKuttaSolver(tableau::AbstractButcherTableau, stepsize::Union{AbstractStepSize, Real}) = ExplicitExponentialRungeKuttaSolver(tableau, stepsize, nothing)
@doc (@doc ExplicitExponentialRungeKuttaSolver) EXRK(args...; kwargs...) = ExplicitExponentialRungeKuttaSolver(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

(solver::ExplicitExponentialRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) = solve!(solution, problem, solver)
(solver::ExplicitExponentialRungeKuttaSolver)(problem::AbstractInitialValueProblem) = solve(problem, solver)
