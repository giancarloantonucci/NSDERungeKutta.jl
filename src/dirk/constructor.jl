# NSDERungeKutta/src/dirk/constructor.jl

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
- `stepsize :: AbstractStepSize`
- `newton :: AbstractNewtonParameters`
- `adaptive :: AbstractAdaptiveParameters` : step-size control for an embedded
  tableau. Note the current limit: a Newton failure inside an adaptive
  implicit step is thrown as a [`NewtonFailure`](@ref) before the controller
  can reject the step. Automatic retry with a shorter step is not implemented
  for implicit solvers in this release; the named DIRK/IRK/IERK solvers are
  fixed-step and are unaffected.

# Methods

    (solver::DiagonallyImplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    (solver::DiagonallyImplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem) :: RungeKuttaSolution

returns the `solution` of a `problem` using `solver`.
"""
struct DiagonallyImplicitRungeKuttaSolver{
    tableau_T<:AbstractButcherTableau,
    stepsize_T<:AbstractStepSize,
    newton_T<:AbstractNewtonParameters,
    adaptive_T<:Union{AbstractAdaptiveParameters,Nothing},
} <: AbstractRungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    newton::newton_T
    adaptive::adaptive_T
    function DiagonallyImplicitRungeKuttaSolver(
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

function DiagonallyImplicitRungeKuttaSolver(
    tableau::AbstractButcherTableau,
    h::Real,
    newton::AbstractNewtonParameters,
    adaptive::Union{AbstractAdaptiveParameters,Nothing},
)
    return DiagonallyImplicitRungeKuttaSolver(
        tableau, StepSize(h; save_stepsizes=false), newton, adaptive
    )
end
function DiagonallyImplicitRungeKuttaSolver(
    tableau::AbstractButcherTableau,
    stepsize::Union{AbstractStepSize,Real},
    newton::AbstractNewtonParameters,
)
    return DiagonallyImplicitRungeKuttaSolver(tableau, stepsize, newton, nothing)
end
@doc (@doc DiagonallyImplicitRungeKuttaSolver) DIRK(args...; kwargs...) =
    DiagonallyImplicitRungeKuttaSolver(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

function (solver::DiagonallyImplicitRungeKuttaSolver)(
    solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem; kwargs...
)
    return solve!(solution, problem, solver; kwargs...)
end
function (solver::DiagonallyImplicitRungeKuttaSolver)(
    problem::AbstractInitialValueProblem; kwargs...
)
    return solve(problem, solver; kwargs...)
end

"Build changed step-size settings while retaining the method and tolerances."
function NSDEBase.withstepsize(solver::DiagonallyImplicitRungeKuttaSolver, h::Real)
    return DiagonallyImplicitRungeKuttaSolver(
        solver.tableau,
        StepSize(h; save_stepsizes=solver.stepsize.save_stepsizes),
        solver.newton,
        solver.adaptive,
    )
end
