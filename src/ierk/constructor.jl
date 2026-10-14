# NSDERungeKutta/src/ierk/constructor.jl

"""
    ImplicitExplicitRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for implicit-explicit solvers.

# Constructors
```julia
ImplicitExplicitRungeKuttaSolver(implicitableau, explicitableau, stepsize, newton[, adaptive])
IERK(args...; kwargs...)
```

# Arguments
- `implicitableau :: AbstractButcherTableau`
- `explicitableau :: AbstractButcherTableau`
- `stepsize :: AbstractStepSize`
- `newton :: AbstractNewtonParameters`
- `adaptive :: AbstractAdaptiveParameters`

# Methods

    (solver::ImplicitExplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    (solver::ImplicitExplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem) :: RungeKuttaSolution

returns the `solution` of a `problem` using `solver`.
"""
struct ImplicitExplicitRungeKuttaSolver{
    implicitableau_T<:AbstractButcherTableau,
    explicitableau_T<:AbstractButcherTableau,
    stepsize_T<:AbstractStepSize,
    newton_T<:AbstractNewtonParameters,
    adaptive_T<:Union{AbstractAdaptiveParameters,Nothing},
} <: AbstractRungeKuttaSolver
    implicitableau::implicitableau_T
    explicitableau::explicitableau_T
    stepsize::stepsize_T
    newton::newton_T
    adaptive::adaptive_T
    function ImplicitExplicitRungeKuttaSolver(
        implicitableau::implicitableau_T,
        explicitableau::explicitableau_T,
        stepsize::stepsize_T,
        newton::newton_T,
        adaptive::adaptive_T,
    ) where {
        implicitableau_T<:AbstractButcherTableau,
        explicitableau_T<:AbstractButcherTableau,
        stepsize_T<:AbstractStepSize,
        newton_T<:AbstractNewtonParameters,
        adaptive_T<:Union{AbstractAdaptiveParameters,Nothing},
    }
        # An IMEX error estimate needs embedded pairs on BOTH tableaus and a
        # two-family controller; neither exists yet, so fail loudly:
        if !(adaptive isa Nothing)
            throw(ArgumentError("Adaptive stepping is not implemented for IMEX solvers in v1."))
        end
        return new{implicitableau_T,explicitableau_T,stepsize_T,newton_T,adaptive_T}(
            implicitableau, explicitableau, stepsize, newton, adaptive
        )
    end
end

function ImplicitExplicitRungeKuttaSolver(
    implicitableau::AbstractButcherTableau,
    explicitableau::AbstractButcherTableau,
    h::Real,
    newton::AbstractNewtonParameters,
    adaptive::Union{AbstractAdaptiveParameters,Nothing},
)
    return ImplicitExplicitRungeKuttaSolver(
        implicitableau, explicitableau, StepSize(h; save_stepsizes=false), newton, adaptive
    )
end
function ImplicitExplicitRungeKuttaSolver(
    implicitableau::AbstractButcherTableau,
    explicitableau::AbstractButcherTableau,
    stepsize::Union{AbstractStepSize,Real},
    newton::AbstractNewtonParameters,
)
    return ImplicitExplicitRungeKuttaSolver(
        implicitableau, explicitableau, stepsize, newton, nothing
    )
end
@doc (@doc ImplicitExplicitRungeKuttaSolver) IERK(args...; kwargs...) =
    ImplicitExplicitRungeKuttaSolver(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

function (solver::ImplicitExplicitRungeKuttaSolver)(
    solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem; kwargs...
)
    return solve!(solution, problem, solver; kwargs...)
end
function (solver::ImplicitExplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem; kwargs...)
    return solve(problem, solver; kwargs...)
end

"Build changed step-size settings while retaining the method and tolerances."
function NSDEBase.withstepsize(solver::ImplicitExplicitRungeKuttaSolver, h::Real)
    return ImplicitExplicitRungeKuttaSolver(
        solver.implicitableau,
        solver.explicitableau,
        StepSize(h; save_stepsizes=solver.stepsize.save_stepsizes),
        solver.newton,
        solver.adaptive,
    )
end
