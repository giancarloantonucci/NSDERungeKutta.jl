abstract type AbstractRungeKuttaSolver <: AbstractInitialValueSolver end
abstract type AbstractRungeKuttaSolution <: AbstractInitialValueSolution end
abstract type AbstractRungeKuttaParameters <: AbstractInitialValueParameters end
abstract type AbstractRungeKuttaCache <: AbstractInitialValueCache end

abstract type AbstractButcherTableau <: AbstractRungeKuttaParameters end
abstract type AbstractStepSize <: AbstractRungeKuttaParameters end
abstract type AbstractAdaptiveParameters <: AbstractRungeKuttaParameters end
abstract type AbstractNewtonParameters <: AbstractRungeKuttaParameters end
