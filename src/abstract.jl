"An abstract type for Runge-Kutta solvers of [`NSDEBase.AbstractInitialValueProblem`](@extref)s."
abstract type AbstractRungeKuttaSolver <: AbstractInitialValueSolver end

"An abstract type for Runge-Kutta solutions of [`NSDEBase.AbstractInitialValueProblem`](@extref)s."
abstract type AbstractRungeKuttaSolution <: AbstractInitialValueSolution end

"An abstract type for Runge-Kutta solver parameters."
abstract type AbstractRungeKuttaParameters <: AbstractInitialValueParameters end

"An abstract type for caching intermediate computations in Runge-Kutta solvers."
abstract type AbstractRungeKuttaCache <: AbstractInitialValueCache end # any cache or intermediate storage

"An abstract type for Butcher tableaus."
abstract type AbstractButcherTableau <: AbstractRungeKuttaParameters end

"An abstract type for time-step sizes, useful for embedded solvers."
abstract type AbstractStepSize <: AbstractRungeKuttaParameters end

"An abstract type for parameters of embedded Runge-Kutta solvers."
abstract type AbstractAdaptiveParameters <: AbstractRungeKuttaParameters end

"An abstract type for parameters of Newton steps in implicit Runge-Kutta solvers."
abstract type AbstractNewtonParameters <: AbstractRungeKuttaParameters end
