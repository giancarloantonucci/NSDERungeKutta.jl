"An abstract type for Runge-Kutta solvers."
abstract type AbstractRungeKuttaSolver <: AbstractInitialValueSolver end

"An abstract type for the solution of an [`AbstractInitialValueProblem`](@ref) obtained with an [`AbstractRungeKuttaSolver`](@ref)."
abstract type AbstractRungeKuttaSolution <: AbstractInitialValueSolution end

"An abstract type for the cache of an [`AbstractRungeKuttaSolver`](@ref)."
abstract type AbstractRungeKuttaCache end

"An abstract type for the parameters of an [`AbstractRungeKuttaSolver`](@ref)."
abstract type AbstractRungeKuttaParameters end

"An abstract type for the Butcher tableau of an [`AbstractRungeKuttaSolver`](@ref)."
abstract type AbstractButcherTableau <: AbstractRungeKuttaParameters end

"An abstract type for the parameters of an adaptive [`AbstractRungeKuttaSolver`](@ref)."
abstract type AbstractAdaptiveParameters <: AbstractRungeKuttaParameters end
