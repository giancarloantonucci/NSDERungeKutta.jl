module RungeKutta

export ButcherTableau
export AdaptiveParameters
export RungeKuttaSolver, RungeKuttaSolution
export ExplicitRungeKuttaSolver, ERK
export ImplicitRungeKuttaSolver, IRK
export stability_function, ℛ

# Fixed-Step Explicit Methods
export Euler, ExplicitEuler
export Midpoint, ExplicitMidpoint
export Heun2
export Ralston2
export Heun3
export Kutta3
export Ralston3
export SSPRK3
export RK4
export Rule38

# Adaptive Explicit Methods
export Fehlberg45, F45
export DormandPrince54, DP54
export Verner65, V65

# Fixed-Step Implicit Methods
export BackwardEuler, ImplicitEuler
export ImplicitMidpoint
export CrankNicolson
export SDIRK3
export HammerHollingsworth, HH4
export LobattoIIIA4
export RadauIIA5

using Reexport
@reexport using NSDEBase
using ArrowMacros
using LinearAlgebra
using BlockArrays
using RecipesBase

zero!(v::AbstractVector) = fill!(v, zero(eltype(v)))
abstract type RungeKuttaSolver <: InitialValueSolver end

include("tableau.jl")
include("adaptive.jl")
include("explicit.jl")
include("implicit.jl")
include("cache.jl")
include("solution.jl")
include("step.jl")
include("solve.jl")
include("stability.jl")
include("recipes.jl")

end
