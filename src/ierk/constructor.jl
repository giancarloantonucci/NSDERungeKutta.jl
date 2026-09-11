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
struct ImplicitExplicitRungeKuttaSolver{implicitableau_T<:AbstractButcherTableau, explicitableau_T<:AbstractButcherTableau, stepsize_T<:AbstractStepSize, newton_T<:AbstractNewtonParameters, adaptive_T<:Union{AbstractAdaptiveParameters,Nothing}} <: AbstractRungeKuttaSolver
    implicitableau :: implicitableau_T
    explicitableau :: explicitableau_T
    stepsize :: stepsize_T
    newton :: newton_T
    adaptive :: adaptive_T
    function ImplicitExplicitRungeKuttaSolver(implicitableau::implicitableau_T, explicitableau::explicitableau_T, stepsize::stepsize_T, newton::newton_T, adaptive::adaptive_T) where {implicitableau_T<:AbstractButcherTableau, explicitableau_T<:AbstractButcherTableau, stepsize_T<:AbstractStepSize, newton_T<:AbstractNewtonParameters, adaptive_T<:Union{AbstractAdaptiveParameters,Nothing}}
        # An IMEX error estimate needs embedded pairs on BOTH tableaus and a
        # two-family controller; neither exists yet, so fail loudly:
        if !(adaptive isa Nothing)
            throw(ArgumentError("Adaptive stepping is not implemented for IMEX solvers in v1."))
        end
        return new{implicitableau_T, explicitableau_T, stepsize_T, newton_T, adaptive_T}(implicitableau, explicitableau, stepsize, newton, adaptive)
    end
end

ImplicitExplicitRungeKuttaSolver(implicitableau::AbstractButcherTableau, explicitableau::AbstractButcherTableau, h::Real, newton::AbstractNewtonParameters, adaptive::Union{AbstractAdaptiveParameters,Nothing}) = ImplicitExplicitRungeKuttaSolver(implicitableau, explicitableau, StepSize(h; save_stepsizes=false), newton, adaptive)
ImplicitExplicitRungeKuttaSolver(implicitableau::AbstractButcherTableau, explicitableau::AbstractButcherTableau, stepsize::Union{AbstractStepSize,Real}, newton::AbstractNewtonParameters) = ImplicitExplicitRungeKuttaSolver(implicitableau, explicitableau, stepsize, newton, nothing)
@doc (@doc ImplicitExplicitRungeKuttaSolver) IERK(args...; kwargs...) = ImplicitExplicitRungeKuttaSolver(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

(solver::ImplicitExplicitRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem; kwargs...) = solve!(solution, problem, solver; kwargs...)
(solver::ImplicitExplicitRungeKuttaSolver)(problem::AbstractInitialValueProblem; kwargs...) = solve(problem, solver; kwargs...)
