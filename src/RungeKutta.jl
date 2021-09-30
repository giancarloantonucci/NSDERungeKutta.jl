module RungeKutta

export ButcherTableau, ℛ
export RungeKuttaSolver, RungeKuttaSolution
export ExplicitRungeKuttaSolver, ERK
export ImplicitRungeKuttaSolver, IRK
export ExplicitExponentialRungeKuttaSolver, EERK
export extract

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

export HeunEuler
export Fehlberg45, F45
export DormandPrince54, DP54
export Verner65, V65

export BackwardEuler, ImplicitEuler
export ImplicitMidpoint
export CrankNicolson
export SDIRK3
export GaussLegendre4, GL4
export GaussLegendre6, GL6
export LobattoIIIA4
export LobattoIIIB2
export LobattoIIIB4
export LobattoIIIC2
export LobattoIIIC4
export RadauIA3
export RadauIA5
export RadauIIA3
export RadauIIA5

export ExponentialRK4, ERK4

using Reexport
using ArrowMacros
using LinearAlgebra
@reexport using NSDEBase
using RecipesBase

"An abstract type for Runge-Kutta solvers."
abstract type RungeKuttaSolver <: InitialValueSolver end

"An abstract type for the cache of a [`RungeKuttaSolver`](@ref)."
abstract type RungeKuttaCache end

include("vector.jl")
include("solution.jl")
include("adaptive.jl")
include("tableau.jl")
include("stepsize.jl")
include("explicit.jl")
include("newton.jl")
include("implicit.jl")
include("exponential.jl")
include("solve.jl")
include("utils.jl")
include("plot.jl")

end
