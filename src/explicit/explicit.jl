"""
    ExplicitRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for an explicit [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
ExplicitRungeKuttaSolver(tableau, stepsize[, adaptive])
ERK(args...; kwargs...)
```

## Arguments
- `tableau :: ButcherTableau`
- `stepsize :: StepSize`
- `adaptive :: AdaptiveParameters` : embedded method's parameters.

# Methods

    (solver::ExplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    
returns the `solution` of a `problem` using `solver`.
    
    (solver::ExplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem; savestages::Bool=false) :: RungeKuttaSolution

returns the [`RungeKuttaSolution`](@ref) of a `problem` using `solver`; `savestages` is a flag to save all stages into `k`.
"""
mutable struct ExplicitRungeKuttaSolver{tableau_T<:ButcherTableau, stepsize_T<:StepSize, adaptive_T<:Union{AdaptiveParameters, Nothing}} <: AbstractRungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    adaptive::adaptive_T
end

ExplicitRungeKuttaSolver(tableau, h::Real, adaptive) = ExplicitRungeKuttaSolver(tableau, StepSize(h), adaptive)
ExplicitRungeKuttaSolver(tableau, stepsize) = ExplicitRungeKuttaSolver(tableau, stepsize, nothing)
@doc (@doc ExplicitRungeKuttaSolver) ERK(args...; kwargs...) = ExplicitRungeKuttaSolver(args...; kwargs...)

#####
##### Methods
#####

(solver::ExplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) = solve!(solution, problem, solver)
(solver::ExplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem; savestages::Bool=false) = solve(problem, solver; savestages=savestages)
