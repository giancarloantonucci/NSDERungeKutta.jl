"""
    ExplicitExponentialRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for an exponential explicit [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
ExplicitExponentialRungeKuttaSolver(tableau, stepsize[, adaptive])
ExRK(args...; kwargs...)
```

# Arguments
- `tableau :: AbstractButcherTableau`
- `stepsize :: AbstractStepSize`
- `adaptive :: AbstractAdaptiveParameters`

# Methods

    (solver::ExplicitExponentialRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    (solver::ExplicitExponentialRungeKuttaSolver)(problem::AbstractInitialValueProblem) :: RungeKuttaSolution

returns the `solution` of a `problem` using `solver`.
"""
struct ExplicitExponentialRungeKuttaSolver{tableau_T<:AbstractButcherTableau, stepsize_T<:AbstractStepSize, adaptive_T<:Union{AbstractAdaptiveParameters,Nothing}} <: AbstractRungeKuttaSolver
    tableau :: tableau_T
    stepsize :: stepsize_T
    adaptive :: adaptive_T
end

ExplicitExponentialRungeKuttaSolver(tableau::AbstractButcherTableau, h::Real, adaptive::Union{AbstractAdaptiveParameters,Nothing}) = ExplicitExponentialRungeKuttaSolver(tableau, StepSize(h; save_stepsizes=false), adaptive)
ExplicitExponentialRungeKuttaSolver(tableau::AbstractButcherTableau, stepsize::Union{AbstractStepSize,Real}) = ExplicitExponentialRungeKuttaSolver(tableau, stepsize, nothing)
@doc (@doc ExplicitExponentialRungeKuttaSolver) ExRK(args...; kwargs...) = ExplicitExponentialRungeKuttaSolver(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

(solver::ExplicitExponentialRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem; kwargs...) = solve!(solution, problem, solver; kwargs...)
(solver::ExplicitExponentialRungeKuttaSolver)(problem::AbstractInitialValueProblem; kwargs...) = solve(problem, solver; kwargs...)
