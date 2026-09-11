# NSDERungeKutta/src/exprk/constructor.jl

"""
    ExponentialRungeKuttaSolver <: AbstractRungeKuttaSolver

A composite type for one-step exponential Runge-Kutta solvers of semilinear
problems `u' = Lu + g(t) + fₙₛ(u, t)`, supplied as a
`NSDEBase.SplitRightHandSide` whose stiff part is a
`NSDEBase.LinearRightHandSide`, or as a plain
`LinearRightHandSide` (on which every scheme propagates the linear flow
exactly, up to the Padé accuracy of the precomputed exponential).

# Constructors
```julia
ExponentialRungeKuttaSolver(tableau, stepsize[, adaptive])
EXPRK(args...; kwargs...)
```

# Arguments
- `tableau :: ExponentialTableau`
- `stepsize :: AbstractStepSize`
- `adaptive :: Nothing` : adaptive stepping is STRUCTURALLY rejected: the
  φ-coefficient operators are evaluated once at `z = h⋅L` when the cache is
  built, which is the entire efficiency argument for these methods at fixed
  `h`; recomputing them on every rejected step would defeat it. Passing
  adaptive parameters throws at construction, mirroring the IMEX solvers.

# Methods

    (solver::ExponentialRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem) :: RungeKuttaSolution
    (solver::ExponentialRungeKuttaSolver)(problem::AbstractInitialValueProblem) :: RungeKuttaSolution

returns the `solution` of a `problem` using `solver`.
"""
struct ExponentialRungeKuttaSolver{tableau_T<:ExponentialTableau, stepsize_T<:AbstractStepSize, adaptive_T<:Union{AbstractAdaptiveParameters,Nothing}} <: AbstractRungeKuttaSolver
    tableau :: tableau_T
    stepsize :: stepsize_T
    adaptive :: adaptive_T
    function ExponentialRungeKuttaSolver(tableau::tableau_T, stepsize::stepsize_T, adaptive::adaptive_T) where {tableau_T<:ExponentialTableau, stepsize_T<:AbstractStepSize, adaptive_T<:Union{AbstractAdaptiveParameters,Nothing}}
        if !(adaptive isa Nothing)
            throw(ArgumentError("Adaptive stepping is not available for exponential solvers: the φ-coefficients are precomputed at z = h⋅L (see the docstring)."))
        end
        return new{tableau_T, stepsize_T, adaptive_T}(tableau, stepsize, adaptive)
    end
end

ExponentialRungeKuttaSolver(tableau::ExponentialTableau, h::Real, adaptive::Union{AbstractAdaptiveParameters,Nothing}) = ExponentialRungeKuttaSolver(tableau, StepSize(h; save_stepsizes=false), adaptive)
ExponentialRungeKuttaSolver(tableau::ExponentialTableau, stepsize::Union{AbstractStepSize,Real}) = ExponentialRungeKuttaSolver(tableau, stepsize, nothing)
@doc (@doc ExponentialRungeKuttaSolver) EXPRK(args...; kwargs...) = ExponentialRungeKuttaSolver(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

(solver::ExponentialRungeKuttaSolver)(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem; kwargs...) = solve!(solution, problem, solver; kwargs...)
(solver::ExponentialRungeKuttaSolver)(problem::AbstractInitialValueProblem; kwargs...) = solve(problem, solver; kwargs...)
