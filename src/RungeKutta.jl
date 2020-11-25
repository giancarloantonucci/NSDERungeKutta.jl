module RungeKutta

export ButcherTableau
export stability_function, ℛ
export AdaptiveParameters
export RungeKuttaSolver, RungeKuttaSolution
export ExplicitRungeKuttaSolver, ERK
export Euler, ExplicitEuler, Midpoint, ExplicitMidpoint, Heun2, Ralston2, Heun3, Kutta3, Ralston3, RK4
export RKF45, RungeKuttaFehlberg, DOPRI54, DormandPrince
export ImplicitRungeKuttaSolver, IRK
export BackwardEuler, ImplicitEuler, ImplicitMidpoint, Trapezoidal, CrankNicolson
export AdditiveRungeKuttaSolver, ARK
export IMEXEuler, AdditiveEuler

using ArrowMacros
using Reexport
@reexport using NSDEBase
using LinearAlgebra
using BlockArrays
using RecipesBase

const ∅ = nothing
zero!(v::AbstractVector) = fill!(v, zero(eltype(v)))
abstract type RungeKuttaSolver <: InitialValueSolver end

include("tableau.jl")
include("stability.jl")
include("adaptive.jl")
include("explicit.jl")
include("implicit.jl")
include("additive.jl")
include("solution.jl")
include("cache.jl")
include("step.jl")
include("solve.jl")
include("recipes.jl")

end
