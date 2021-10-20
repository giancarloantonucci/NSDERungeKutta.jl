"""
    ImplicitRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for an implicit [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
ImplicitRungeKuttaSolver(tableau, h, newton[, adaptive])
IRK(args...; kwargs...)
```

# Arguments
- `tableau  :: ButcherTableau`     : Butcher tableau.
- `stepsize :: StepSize`           : step-size.
- `newton   :: NewtonParameters`   : simplified Newton's parameters.
- `adaptive :: AdaptiveParameters` : embedded method's parameters.

# Functions
- [`show`   ](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
mutable struct ImplicitRungeKuttaSolver{tableau_T, stepsize_T, newton_T, adaptive_T} <: AbstractRungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    newton::newton_T
    adaptive::adaptive_T
end

ImplicitRungeKuttaSolver(tableau, h::Real, newton, adaptive) = ImplicitRungeKuttaSolver(tableau, StepSize(h), newton, adaptive)
ImplicitRungeKuttaSolver(tableau, stepsize, newton) = ImplicitRungeKuttaSolver(tableau, stepsize, newton, nothing)
@doc (@doc ImplicitRungeKuttaSolver) IRK(args...; kwargs...) = ImplicitRungeKuttaSolver(args...; kwargs...)

############################################################################################
#                                         PRINTING                                         #
############################################################################################

"""
    show(io::IO, solver::ImplicitRungeKuttaSolver)

prints a full description of `solver` and its contents to a stream `io`.
"""
Base.show(io::IO, solver::ImplicitRungeKuttaSolver) = NSDEBase._show(io, solver)

"""
    summary(io::IO, solver::ImplicitRungeKuttaSolver)

prints a brief description of `solver` to a stream `io`.
"""
Base.summary(io::IO, solver::ImplicitRungeKuttaSolver) = NSDEBase._summary(io, solver)

############################################################################################
#                                          METHODS                                         #
############################################################################################

function (solver::ImplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem; save_stages::Bool = false)
    solve(problem, solver; save_stages=save_stages)
end
