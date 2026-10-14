# NSDERungeKutta/src/erk/constructor.jl

"""
    ExplicitRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for explicit solvers.

# Constructors
```julia
ExplicitRungeKuttaSolver(tableau, stepsize[, adaptive])
ERK(args...; kwargs...)
```

# Arguments
- `tableau :: AbstractButcherTableau`
- `stepsize :: AbstractStepSize`
- `adaptive :: AbstractAdaptiveParameters`

# Methods

    (solver::ExplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    (solver::ExplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem) :: RungeKuttaSolution

returns the `solution` of a `problem` using `solver`.
"""
struct ExplicitRungeKuttaSolver{
    tableau_T<:AbstractButcherTableau,
    stepsize_T<:AbstractStepSize,
    adaptive_T<:Union{AbstractAdaptiveParameters,Nothing},
} <: AbstractRungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    adaptive::adaptive_T
    function ExplicitRungeKuttaSolver(
        tableau::tableau_T, stepsize::stepsize_T, adaptive::adaptive_T
    ) where {
        tableau_T<:AbstractButcherTableau,
        stepsize_T<:AbstractStepSize,
        adaptive_T<:Union{AbstractAdaptiveParameters,Nothing},
    }
        check_adaptive(tableau, adaptive) # adaptivity needs an embedded pair; fail at construction
        return new{tableau_T,stepsize_T,adaptive_T}(tableau, stepsize, adaptive)
    end
end

function ExplicitRungeKuttaSolver(
    tableau::AbstractButcherTableau, h::Real, adaptive::Union{AbstractAdaptiveParameters,Nothing}
)
    return ExplicitRungeKuttaSolver(tableau, StepSize(h; save_stepsizes=false), adaptive)
end
function ExplicitRungeKuttaSolver(
    tableau::AbstractButcherTableau, stepsize::Union{AbstractStepSize,Real}
)
    return ExplicitRungeKuttaSolver(tableau, stepsize, nothing)
end
@doc (@doc ExplicitRungeKuttaSolver) ERK(args...; kwargs...) =
    ExplicitRungeKuttaSolver(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

function (solver::ExplicitRungeKuttaSolver)(
    solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem; kwargs...
)
    return solve!(solution, problem, solver; kwargs...)
end
function (solver::ExplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem; kwargs...)
    return solve(problem, solver; kwargs...)
end

"Build changed step-size settings while retaining the method and tolerances."
function NSDEBase.withstepsize(solver::ExplicitRungeKuttaSolver, h::Real)
    return ExplicitRungeKuttaSolver(
        solver.tableau, StepSize(h; save_stepsizes=solver.stepsize.save_stepsizes), solver.adaptive
    )
end
