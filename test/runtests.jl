using LinearAlgebra
using NSDERungeKutta
using Test
# using Aqua

# Aqua.test_all(NSDEBase)

u0 = ones(3)
tspan = (0.0, 1e-3)
function f!(du, u, t)
    I₃ = Matrix(1.0I, 3, 3)
    du .=  I₃ * u
end
problem = IVP(f!, u0, tspan)

@testset "ERK" begin
    solver = ExplicitEuler(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Heun2(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = ExplicitMidpoint(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Ralston2(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Heun3(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = RK3(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Ralston3(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = SSPRK3(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Ralston4(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = RK4(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Rule38(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Butcher5(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = KuttaNystrom5(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Butcher6(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Butcher7(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = HeunEuler(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = F45(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = DP54(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = V65(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = F78(h = 1e-3)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
end

@testset "DIRK" begin
    solver = ImplicitEuler(h = 1e-3)
    @test solver isa DiagonallyImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = GaussLegendre2(h = 1e-3)
    @test solver isa DiagonallyImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = SDIRK2(h = 1e-3)
    @test solver isa DiagonallyImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = LobattoIII2(h = 1e-3)
    @test solver isa DiagonallyImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = LobattoIIIA2(h = 1e-3)
    @test solver isa DiagonallyImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = SDIRK3(h = 1e-3)
    @test solver isa DiagonallyImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = RadauI3(h = 1e-3)
    @test solver isa DiagonallyImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = RadauII3(h = 1e-3)
    @test solver isa DiagonallyImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = SDIRK4(h = 1e-3)
    @test solver isa DiagonallyImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = LobattoIII4(h = 1e-3)
    @test solver isa DiagonallyImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
end

@testset "IRK" begin
    solver = LobattoIIIC2(h = 1e-3)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = RadauIA3(h = 1e-3)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = RadauIIA3(h = 1e-3)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = GaussLegendre4(h = 1e-3)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = LobattoIIIA4(h = 1e-3)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = LobattoIIIB4(h = 1e-3)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = LobattoIIIC4(h = 1e-3)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = RadauI5(h = 1e-3)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = RadauIA5(h = 1e-3)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = RadauII5(h = 1e-3)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = RadauIIA5(h = 1e-3)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = GaussLegendre6(h = 1e-3)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
end

@testset "ℛ(z)" begin
    @test ℛ(0.0, Euler().tableau) == 1.0
end
