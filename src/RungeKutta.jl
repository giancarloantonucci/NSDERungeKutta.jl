module RungeKutta

export ButcherTableau
export RungeKuttaSolver, RungeKuttaSolution
export ExplicitRungeKuttaSolver, ERK
export ImplicitRungeKuttaSolver, IRK
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

using Reexport
using ArrowMacros
using LinearAlgebra
@reexport using NSDEBase
using RecipesBase

abstract type RungeKuttaSolver <: InitialValueSolver end

Vector{T}(undef, n, d) where T = Vector{T}[Vector{T}(undef, d) for i = 1:n]
Vector{T}(undef, m, n, d) where T = Vector{Vector{T}}[Vector{T}(undef, n, d) for i = 1:m]

zero!(v::AbstractVector) = fill!(v, zero(eltype(v)))
function zero!(v::AbstractVector{<:AbstractVector})
    for i in eachindex(v)
        zero!(v[i])
    end
    return v
end

function norm!(v::AbstractVector{<:AbstractVector})
    r = zero(eltype(v))
    for i in eachindex(v)
        r += norm(v[i])
    end
    return r
end

include("tableau.jl")
include("stepsize.jl")
include("adaptive.jl")
include("explicit.jl")
include("implicit.jl")
include("cache.jl")
include("solution.jl")
include("step.jl")
include("solve.jl")
include("stability.jl")
include("recipes.jl")
# include("show.jl")

end
