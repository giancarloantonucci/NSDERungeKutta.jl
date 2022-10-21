module RungeKutta

using Reexport
using ArrowMacros
using LinearAlgebra
@reexport using NSDEBase
using RecipesBase

include("abstract.jl")
include("utils.jl")
include("solution.jl")
include("solve.jl")

include("tableau.jl")
include("stepsize.jl")
include("newton.jl")
include("adaptive.jl")

include("erk/constructor.jl")
include("erk/cache.jl")
include("erk/step.jl")
include("erk/adaptivestep.jl")
include("erk/solvers.jl")

include("dirk/constructor.jl")
include("dirk/cache.jl")
include("dirk/step.jl")
include("dirk/adaptivestep.jl")
include("dirk/solvers.jl")

include("ierk/constructor.jl")
include("ierk/cache.jl")
include("ierk/step.jl")
include("ierk/adaptivestep.jl")
include("ierk/solvers.jl")

include("exrk/constructor.jl")
include("exrk/cache.jl")
include("exrk/step.jl")
include("exrk/adaptivestep.jl")
include("exrk/solvers.jl")

include("irk/constructor.jl")
include("irk/cache.jl")
include("irk/step.jl")
include("irk/adaptivestep.jl")
include("irk/solvers.jl")

include("stability.jl")
include("plot.jl")

export AbstractRungeKuttaSolver
export AbstractRungeKuttaSolution
export AbstractRungeKuttaParameters

export RungeKuttaSolution
export ButcherTableau

export ExplicitRungeKuttaSolver, ERK
export Euler, ExplicitEuler
export Heun2
export Midpoint, ExplicitMidpoint
export Ralston2
export Heun3
export Kutta3
export Ralston3
export SSPRK3
export Ralston4
export RungeKutta4, RK4
export Rule38
export Butcher5
export KuttaNyström5
export Butcher6
export Butcher7

export HeunEuler
export BogackiShampine
export Fehlberg45
export DormandPrince54
export Verner65
export Fehlberg78

export DiagonallyImplicitRungeKuttaSolver, DIRK
export BackwardEuler, ImplicitEuler
export ImplicitMidpoint, GaußLegendre2, GaussLegendre2
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

export ExplicitExponentialRungeKuttaSolver, EXRK
export ExponentialRK4, ERK4

export ImplicitRungeKuttaSolver, IRK
export LobattoIIIC2
export RadauIA3
export RadauIIA3
export GaußLegendre4, GaussLegendre4
export LobattoIIIA4
export LobattoIIIB4
export LobattoIIIC4
export RadauI5
export RadauIA5
export RadauII5
export RadauIIA5
export GaußLegendre6, GaussLegendre6

export butchertableau
export stepsize
export dimension, extract
export ℛ

end
