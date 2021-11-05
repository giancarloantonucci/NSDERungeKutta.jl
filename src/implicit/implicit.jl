"""
    ImplicitRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for an implicit [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
ImplicitRungeKuttaSolver(tableau, stepsize, newton[, adaptive])
IRK(args...; kwargs...)
FIRK(args...; kwargs...)
```

## Arguments
- `tableau :: ButcherTableau`
- `stepsize :: StepSize`
- `newton :: NewtonParameters` : simplified Newton's parameters.
- `adaptive :: AdaptiveParameters` : embedded method's parameters.

# Methods

    (solver::ImplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    
returns the `solution` of a `problem` using `solver`.
    
    (solver::ImplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem; savestages::Bool=false) :: RungeKuttaSolution

returns the [`RungeKuttaSolution`](@ref) of a `problem` using `solver`; `savestages` is a flag to save all stages into `k`.
"""
mutable struct ImplicitRungeKuttaSolver{tableau_T<:ButcherTableau, stepsize_T<:StepSize, newton_T<:NewtonParameters, adaptive_T<:Union{AdaptiveParameters, Nothing}} <: AbstractRungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    newton::newton_T
    adaptive::adaptive_T
end

ImplicitRungeKuttaSolver(tableau, h::Real, newton, adaptive) = ImplicitRungeKuttaSolver(tableau, StepSize(h), newton, adaptive)
ImplicitRungeKuttaSolver(tableau, stepsize, newton) = ImplicitRungeKuttaSolver(tableau, stepsize, newton, nothing)
@doc (@doc ImplicitRungeKuttaSolver) IRK(args...; kwargs...) = ImplicitRungeKuttaSolver(args...; kwargs...)
@doc (@doc ImplicitRungeKuttaSolver) FIRK(args...; kwargs...) = ImplicitRungeKuttaSolver(args...; kwargs...)

#####
##### Methods
#####

(solver::ImplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) = solve!(solution, problem, solver)
(solver::ImplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem; savestages::Bool=false) = solve(problem, solver; savestages=savestages)
