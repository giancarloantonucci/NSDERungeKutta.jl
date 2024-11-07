"An abstract type for a Runge-Kutta solver."
abstract type AbstractRungeKuttaSolver <: AbstractInitialValueSolver end

"An abstract type for a Runge-Kutta solution."
abstract type AbstractRungeKuttaSolution <: AbstractInitialValueSolution end

"An abstract type for Runge-Kutta solver parameters."
abstract type AbstractRungeKuttaParameters <: AbstractInitialValueParameters end

"An abstract type for the temporary variables in a Runge-Kutta solver."
abstract type AbstractRungeKuttaCache <: AbstractInitialValueCache end # any cache or intermediate storage

"An abstract type for a Butcher tableau."
abstract type AbstractButcherTableau <: AbstractRungeKuttaParameters end

"An abstract type for the time-step size, useful for adaptive solvers."
abstract type AbstractStepSize <: AbstractRungeKuttaParameters end

"An abstract type for the parameters of an adaptive (or embedded) Runge-Kutta solver."
abstract type AbstractAdaptiveParameters <: AbstractRungeKuttaParameters end

"An abstract type for parameters of the Newton steps in an implicit Runge-Kutta solver."
abstract type AbstractNewtonParameters <: AbstractRungeKuttaParameters end
