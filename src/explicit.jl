"""
    ExplicitRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for an explicit [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
ExplicitRungeKuttaSolver(tableau, h[, adaptive])
ERK(args...; kwargs...)
```

# Arguments
- `tableau  :: ButcherTableau`     : Butcher tableau.
- `stepsize :: StepSize`           : step-size.
- `adaptive :: AdaptiveParameters` : embedded method's parameters.

# Functions
- [`show`   ](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
mutable struct ExplicitRungeKuttaSolver{tableau_T, stepsize_T, adaptive_T} <: AbstractRungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    adaptive::adaptive_T
end

ExplicitRungeKuttaSolver(tableau, h::Real, adaptive) = ExplicitRungeKuttaSolver(tableau, StepSize(h), adaptive)
ExplicitRungeKuttaSolver(tableau, stepsize) = ExplicitRungeKuttaSolver(tableau, stepsize, nothing)
@doc (@doc ExplicitRungeKuttaSolver) ERK(args...; kwargs...) = ExplicitRungeKuttaSolver(args...; kwargs...)

############################################################################################
#                                         PRINTING                                         #
############################################################################################

"""
    show(io::IO, solver::ExplicitRungeKuttaSolver)

prints a full description of `solver` and its contents to a stream `io`.
"""
Base.show(io::IO, solver::ExplicitRungeKuttaSolver) = NSDEBase._show(io, solver)

"""
    summary(io::IO, solver::ExplicitRungeKuttaSolver)

prints a brief description of `solver` to a stream `io`.
"""
Base.summary(io::IO, solver::ExplicitRungeKuttaSolver) = NSDEBase._summary(io, solver)

############################################################################################
#                                          METHODS                                         #
############################################################################################

function (solver::ExplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem; save_stages::Bool = false)
    solve(problem, solver; save_stages=save_stages)
end
