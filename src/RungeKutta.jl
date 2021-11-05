module RungeKutta

using ArrowMacros
using LinearAlgebra
using NSDEBase
using RecipesBase

include("abstract.jl")
include("utils.jl")
include("tmp.jl")
include("solution.jl")
include("solve.jl")
include("tableau.jl")
include("stepsize.jl")
include("adaptive.jl")
include("newton.jl")
include("explicit/explicit.jl")
include("explicit/cache.jl")
include("explicit/step.jl")
include("explicit/solvers.jl")
include("implicit/implicit.jl")
include("implicit/cache.jl")
include("implicit/step.jl")
include("implicit/solvers.jl")
include("exponential/exponential.jl")
include("exponential/cache.jl")
include("exponential/step.jl")
include("exponential/solvers.jl")
include("stability.jl")
include("plotrecipes.jl")

export AbstractRungeKuttaSolver
export AbstractRungeKuttaSolution
export AbstractRungeKuttaCache
export AbstractRungeKuttaParameters
export AbstractAdaptiveParameters
export AbstractButcherTableau
export AbstractStepSize
export AbstractNewtonParameters

export ExplicitRungeKuttaSolver, ERK
export ImplicitRungeKuttaSolver, IRK, FIRK
export ExplicitExponentialRungeKuttaSolver, EERK
export RungeKuttaSolution
export ButcherTableau
export StepSize
export AdaptiveParameters
export NewtonParameters

export butchertableau
export stepsize
export adaptivecheck!
export extract
export ℛ

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

end
