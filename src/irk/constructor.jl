"""
    ImplicitRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for implicit solvers.

# Constructors
```julia
ImplicitRungeKuttaSolver(tableau, stepsize, newton[, adaptive])
IRK(args...; kwargs...)
```

# Arguments
- `tableau :: AbstractButcherTableau`
- `stepsize :: Union{AbstractStepSize, Real}`
- `newton :: AbstractNewtonParameters`
- `adaptive :: Union{AbstractAdaptiveParameters, Nothing}`

# Methods

    (solver::ImplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    (solver::ImplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem) :: RungeKuttaSolution

returns the `solution` of a `problem` using `solver`.
"""
struct ImplicitRungeKuttaSolver{tableau_T<:AbstractButcherTableau, stepsize_T<:AbstractStepSize, newton_T<:AbstractNewtonParameters, adaptive_T<:Union{AbstractAdaptiveParameters, Nothing}} <: AbstractRungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    newton::newton_T
    adaptive::adaptive_T
end

ImplicitRungeKuttaSolver(tableau::AbstractButcherTableau, h::Real, newton::AbstractNewtonParameters, adaptive::Union{AbstractAdaptiveParameters, Nothing}) = ImplicitRungeKuttaSolver(tableau, StepSize(h), newton, adaptive)
ImplicitRungeKuttaSolver(tableau::AbstractButcherTableau, stepsize::Union{AbstractStepSize, Real}, newton::AbstractNewtonParameters) = ImplicitRungeKuttaSolver(tableau, stepsize, newton, nothing)
@doc (@doc ImplicitRungeKuttaSolver) IRK(args...; kwargs...) = ImplicitRungeKuttaSolver(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

(solver::ImplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) = solve!(solution, problem, solver)
(solver::ImplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem) = solve(problem, solver)
