using NSDERungeKutta
using LinearAlgebra
using Test
using Aqua
using RecipesBase: RecipesBase
using NSDEBase: NSDEBase

include("support.jl")

@testset "NSDERungeKutta" begin
    @testset "Aqua" begin
        # Same reasoning as NSDEBase: recipetype(::Val{…}) is the established
        # user-plot idiom; whitelist that one function, keep the rest strict.
        Aqua.test_all(
            NSDERungeKutta; piracies=(; treat_as_own=[NSDERungeKutta.RecipesBase.recipetype])
        )
    end

    include("tableaus_and_stability.jl")
    include("convergence.jl")
    include("order_conditions.jl")
    include("reliability.jl")
    include("adaptive.jl")
    include("adaptive_controller.jl")
    include("numeric_types.jl")
    include("solutions.jl")
    include("newton.jl")
    include("endpoints.jl")
    include("allocations.jl")
    include("exponential_methods.jl")
    include("exponential_oracles.jl")
end
