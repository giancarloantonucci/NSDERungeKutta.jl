"""
    ImplicitExplicitRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for implicit-explicit solvers.

# Constructors
```julia
ImplicitExplicitRungeKuttaSolver(implicitableau, explicitableau, stepsize, newton[, adaptive])
IERK(args...; kwargs...)
```

# Arguments
- `implicitableau :: AbstractButcherTableau`
- `explicitableau :: AbstractButcherTableau`
- `stepsize :: Union{AbstractStepSize, Real}`
- `newton :: AbstractNewtonParameters`
- `adaptive :: Union{AbstractAdaptiveParameters, Nothing}`

# Methods

    (solver::ImplicitExplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    (solver::ImplicitExplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem) :: RungeKuttaSolution

returns the `solution` of a `problem` using `solver`.
"""
struct ImplicitExplicitRungeKuttaSolver{implicitableau_T<:AbstractButcherTableau, explicitableau_T<:AbstractButcherTableau, stepsize_T<:AbstractStepSize, newton_T<:AbstractNewtonParameters, adaptive_T<:Union{AbstractAdaptiveParameters, Nothing}} <: AbstractRungeKuttaSolver
    implicitableau::implicitableau_T
    explicitableau::explicitableau_T
    stepsize::stepsize_T
    newton::newton_T
    adaptive::adaptive_T
end

ImplicitExplicitRungeKuttaSolver(implicitableau::AbstractButcherTableau, explicitableau::AbstractButcherTableau, h::Real, newton::AbstractNewtonParameters, adaptive::Union{AbstractAdaptiveParameters, Nothing}) = ImplicitExplicitRungeKuttaSolver(implicitableau, explicitableau, StepSize(h), newton, adaptive)
ImplicitExplicitRungeKuttaSolver(implicitableau::AbstractButcherTableau, explicitableau::AbstractButcherTableau, stepsize::Union{AbstractStepSize, Real}, newton::AbstractNewtonParameters) = ImplicitExplicitRungeKuttaSolver(implicitableau, explicitableau, stepsize, newton, nothing)
@doc (@doc ExplicitRungeKuttaSolver) IERK(args...; kwargs...) = ImplicitExplicitRungeKuttaSolver(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

(solver::ImplicitExplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) = solve!(solution, problem, solver)
(solver::ImplicitExplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem) = solve(problem, solver)
