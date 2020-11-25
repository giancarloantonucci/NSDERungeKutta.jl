using RungeKutta
using Test

@testset "Explicit" begin
    problem = Lorenz()
    solver = ExplicitEuler(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solver = ExplicitMidpoint(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solver = Heun2(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solver = Ralston2(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solver = Heun3(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solver = Kutta3(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solver = Ralston3(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solver = RK4(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solver = RKF45(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solver = DOPRI54(h = 1e-2)
    @test solver isa ExplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
end

@testset "Implicit" begin
    problem = Lorenz()
    solver = ImplicitEuler(h = 1e-2)
    @test solver isa ImplicitRungeKuttaSolver
    solver = ImplicitMidpoint(h = 1e-2)
    @test solver isa ImplicitRungeKuttaSolver
    solver = CrankNicolson(h = 1e-2)
    @test solver isa ImplicitRungeKuttaSolver
    solution = solve(problem, solver)
    @test solution isa RungeKuttaSolution
end

@testset "Stability" begin
    @test ℛ(0.0, Euler().tableau) == 1.0
end
