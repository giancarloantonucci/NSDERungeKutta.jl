"""
    DiagonallyImplicitRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for diagonally-implicit solvers.

# Constructors
```julia
DiagonallyImplicitRungeKuttaSolver(tableau, stepsize, newton[, adaptive])
DIRK(args...; kwargs...)
```

# Arguments
- `tableau :: AbstractButcherTableau`
- `stepsize :: Union{AbstractStepSize, Real}`
- `newton :: AbstractNewtonParameters`
- `adaptive :: Union{AbstractAdaptiveParameters, Nothing}`

# Methods

    (solver::DiagonallyImplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    (solver::DiagonallyImplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem) :: RungeKuttaSolution

returns the `solution` of a `problem` using `solver`.
"""
struct DiagonallyImplicitRungeKuttaSolver{tableau_T<:AbstractButcherTableau, stepsize_T<:AbstractStepSize, newton_T<:AbstractNewtonParameters, adaptive_T<:Union{AbstractAdaptiveParameters, Nothing}} <: AbstractRungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    newton::newton_T
    adaptive::adaptive_T
end

DiagonallyImplicitRungeKuttaSolver(tableau::AbstractButcherTableau, h::Real, newton::AbstractNewtonParameters, adaptive::Union{AbstractAdaptiveParameters, Nothing}) = DiagonallyImplicitRungeKuttaSolver(tableau, StepSize(h), newton, adaptive)
DiagonallyImplicitRungeKuttaSolver(tableau::AbstractButcherTableau, stepsize::Union{AbstractStepSize, Real}, newton::AbstractNewtonParameters) = DiagonallyImplicitRungeKuttaSolver(tableau, stepsize, newton, nothing)
@doc (@doc DiagonallyImplicitRungeKuttaSolver) DIRK(args...; kwargs...) = DiagonallyImplicitRungeKuttaSolver(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

(solver::DiagonallyImplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) = solve!(solution, problem, solver)
(solver::DiagonallyImplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem) = solve(problem, solver)
