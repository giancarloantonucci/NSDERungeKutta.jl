using LinearAlgebra
using RungeKutta
using Test

u0 = [1.0, 1.0, 1.0]
tspan = (0.0, 1e-2)
function f!(du, u, t)
    du .= Matrix(1.0I, 3, 3) * u
end
problem = IVP(f!, u0, tspan)

@testset "Explicit" begin
    solver = ExplicitEuler(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = ExplicitMidpoint(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Heun2(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Ralston2(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Heun3(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Kutta3(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Ralston3(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = RK4(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = Rule38(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = F45(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = DP54(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = V65(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
end

@testset "Implicit" begin
    solver = ImplicitEuler(h = 1e-2)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = ImplicitMidpoint(h = 1e-2)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = CrankNicolson(h = 1e-2)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = SDIRK3(h = 1e-2)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = HH4(h = 1e-2)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = LobattoIIIA4(h = 1e-2)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
    solver = RadauIIA5(h = 1e-2)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
end

@testset "Stability" begin
    @test ℛ(0.0, Euler().tableau) == 1.0
end
