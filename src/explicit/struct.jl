"""
    ExplicitRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for an explicit [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
ExplicitRungeKuttaSolver(tableau, stepsize[, adaptive])
ERK(args...; kwargs...)
```

# Arguments
- `tableau :: ButcherTableau` : Butcher tableau.
- `stepsize :: StepSize` : step-size.
- `adaptive :: AdaptiveParameters` : embedded method's parameters.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
mutable struct ExplicitRungeKuttaSolver{tableau_T, stepsize_T, adaptive_T} <: AbstractRungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    adaptive::adaptive_T
end

function ExplicitRungeKuttaSolver(tableau, h::Real, adaptive)
    return ExplicitRungeKuttaSolver(tableau, StepSize(h), adaptive)
end

function ExplicitRungeKuttaSolver(tableau, stepsize)
    return ExplicitRungeKuttaSolver(tableau, stepsize, nothing)
end

@doc (@doc ExplicitRungeKuttaSolver) function ERK(args...; kwargs...)
    return ExplicitRungeKuttaSolver(args...; kwargs...)
end

#####
##### Methods
#####

# function (solver::ExplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem; savestages::Bool = false)
#     return solve!(solution, problem, solver; savestages=savestages)
# end

# function (solver::ExplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem; savestages::Bool = false)
#     return solve(problem, solver; savestages=savestages)
# end
