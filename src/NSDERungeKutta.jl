# NSDERungeKutta/src/NSDERungeKutta.jl

module NSDERungeKutta

using Reexport
using ArrowMacros
using LinearAlgebra
import SparseArrays
@reexport using NSDEBase
using RecipesBase

include("abstract.jl")
include("utils.jl")
include("solution.jl")
include("solve.jl")

include("tableau.jl")
include("tableau_show.jl")
include("stepsize.jl")
include("newton.jl")
include("adaptive.jl")

include("erk/constructor.jl")
include("erk/cache.jl")
include("erk/step.jl")
include("erk/solvers.jl")

include("dirk/constructor.jl")
include("dirk/cache.jl")
include("dirk/step.jl")
include("dirk/solvers.jl")

include("ierk/constructor.jl")
include("ierk/cache.jl")
include("ierk/step.jl")
include("ierk/solvers.jl")

include("irk/constructor.jl")
include("irk/cache.jl")
include("irk/step.jl")
include("irk/solvers.jl")

include("phi.jl")
include("exprk/tableau.jl")
include("exprk/constructor.jl")
include("exprk/cache.jl")
include("exprk/step.jl")
include("exprk/solvers.jl")

include("stability.jl")
include("plots_recipes.jl")

export AbstractRungeKuttaSolver
export AbstractRungeKuttaSolution
export AbstractRungeKuttaParameters

export RungeKuttaSolution
export ButcherTableau
export NewtonParameters, NewtonFailure

export ExplicitRungeKuttaSolver, ERK
export Euler, ExplicitEuler
export Heun2
export Midpoint, ExplicitMidpoint
export Ralston2
export Heun3
export RungeKutta3, RK3
export Ralston3
export SSPRK3
export Ralston4
export RungeKutta4, RK4
export Rule38
export Butcher5
export KuttaNystrom5
export Butcher6
export Butcher7

export HeunEuler
export BogackiShampine
export Fehlberg45, F45
export DormandPrince54, DP54
export Verner65, V65
export Fehlberg78, F78

export DiagonallyImplicitRungeKuttaSolver, DIRK
export BackwardEuler, ImplicitEuler
export ImplicitMidpoint, GaussLegendre2
export SDIRK2
export LobattoIII2
export CrankNicolson, LobattoIIIA2
export SDIRK3
export RadauI3
export RadauII3
export SDIRK4
export LobattoIII4

export ImplicitExplicitRungeKuttaSolver, IERK
export IMEXEuler, IMEXSSP1_111
export IMEXSSP2_222
export IMEXSSP2_322
export IMEXSSP2_332
export IMEXSSP3_332

export ExponentialRungeKuttaSolver, EXPRK
export ExponentialTableau
export phifunctions, expphifunctions
export LawsonEuler, Lawson4
export NorsettEuler, ETDEuler, ExponentialEuler
export ETD2RK, ETD3RK, ETD4RK, ETDRK4
export Krogstad
export HochbruckOstermann4, HochOst4

export ImplicitRungeKuttaSolver, IRK
export LobattoIIIC2
export RadauIA3
export RadauIIA3
export GaussLegendre4
export LobattoIIIA4
export LobattoIIIB4
export LobattoIIIC4
export RadauI5
export RadauIA5
export RadauII5
export RadauIIA5
export GaussLegendre6

export butchertableau
export stepsize
export numtimesteps, numvariables, extract
export stability_function

end
