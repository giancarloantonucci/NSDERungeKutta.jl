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
- `stepsize :: AbstractStepSize`
- `newton :: AbstractNewtonParameters`
- `adaptive :: AbstractAdaptiveParameters`

# Methods

    (solver::ImplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    (solver::ImplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem) :: RungeKuttaSolution

returns the `solution` of a `problem` using `solver`.
"""
struct ImplicitRungeKuttaSolver{
    tableau_T<:AbstractButcherTableau,
    stepsize_T<:AbstractStepSize,
    newton_T<:AbstractNewtonParameters,
    adaptive_T<:Union{AbstractAdaptiveParameters,Nothing},
} <: AbstractRungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    newton::newton_T
    adaptive::adaptive_T
    function ImplicitRungeKuttaSolver(
        tableau::tableau_T, stepsize::stepsize_T, newton::newton_T, adaptive::adaptive_T
    ) where {
        tableau_T<:AbstractButcherTableau,
        stepsize_T<:AbstractStepSize,
        newton_T<:AbstractNewtonParameters,
        adaptive_T<:Union{AbstractAdaptiveParameters,Nothing},
    }
        check_adaptive(tableau, adaptive) # adaptivity needs an embedded pair; fail at construction
        return new{tableau_T,stepsize_T,newton_T,adaptive_T}(tableau, stepsize, newton, adaptive)
    end
end

function ImplicitRungeKuttaSolver(
    tableau::AbstractButcherTableau,
    h::Real,
    newton::AbstractNewtonParameters,
    adaptive::Union{AbstractAdaptiveParameters,Nothing},
)
    return ImplicitRungeKuttaSolver(tableau, StepSize(h; save_stepsizes=false), newton, adaptive)
end
function ImplicitRungeKuttaSolver(
    tableau::AbstractButcherTableau,
    stepsize::Union{AbstractStepSize,Real},
    newton::AbstractNewtonParameters,
)
    return ImplicitRungeKuttaSolver(tableau, stepsize, newton, nothing)
end
@doc (@doc ImplicitRungeKuttaSolver) IRK(args...; kwargs...) =
    ImplicitRungeKuttaSolver(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

function (solver::ImplicitRungeKuttaSolver)(
    solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem; kwargs...
)
    return solve!(solution, problem, solver; kwargs...)
end
function (solver::ImplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem; kwargs...)
    return solve(problem, solver; kwargs...)
end

"Build changed step-size settings while retaining the method and tolerances."
function NSDEBase.withstepsize(solver::ImplicitRungeKuttaSolver, h::Real)
    return ImplicitRungeKuttaSolver(
        solver.tableau,
        StepSize(h; save_stepsizes=solver.stepsize.save_stepsizes),
        solver.newton,
        solver.adaptive,
    )
end
