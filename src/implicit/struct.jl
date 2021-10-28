"""
    ImplicitRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for an implicit [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
ImplicitRungeKuttaSolver(tableau, stepsize, newton[, adaptive])
IRK(args...; kwargs...)
```

# Arguments
- `tableau :: ButcherTableau` : Butcher tableau.
- `stepsize :: StepSize` : step-size.
- `newton :: NewtonParameters` : simplified Newton's parameters.
- `adaptive :: AdaptiveParameters` : embedded method's parameters.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
mutable struct ImplicitRungeKuttaSolver{tableau_T, stepsize_T, newton_T, adaptive_T} <: AbstractRungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    newton::newton_T
    adaptive::adaptive_T
end

function ImplicitRungeKuttaSolver(tableau, h::Real, newton, adaptive)
    return ImplicitRungeKuttaSolver(tableau, StepSize(h), newton, adaptive)
end

function ImplicitRungeKuttaSolver(tableau, stepsize, newton)
    return ImplicitRungeKuttaSolver(tableau, stepsize, newton, nothing)
end

@doc (@doc ImplicitRungeKuttaSolver) function IRK(args...; kwargs...)
    return ImplicitRungeKuttaSolver(args...; kwargs...)
end

#####
##### Methods
#####

# function (solver::ImplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem; savestages::Bool = false)
#     return solve!(solution, problem, solver; savestages=savestages)
# end

# function (solver::ImplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem; savestages::Bool = false)
#     return solve(problem, solver; savestages=savestages)
# end
