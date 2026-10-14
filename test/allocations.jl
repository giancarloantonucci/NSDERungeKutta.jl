@testset "allocation-free ERK step" begin
    L = [0.0 1.0; -4.0 -0.4]
    g!(dg, t) = (dg[1]=sin(t); dg[2]=cos(t); dg)
    f!(du, u, t) = (@. du = -u; du)
    problems = (
        IVP(f!, ones(2), (0.0, 1.0)),                                     # NonlinearRHS
        IVP(LRHS(L, g!), ones(2), (0.0, 1.0)),                            # LRHS with forcing
        IVP(SRHS(LRHS(L, g!), RHS((du, u, t) -> (@. du = sin(u); du))), ones(2), (0.0, 1.0)), # split
    )
    for problem in problems
        @test step_allocs(problem, RungeKutta4(h=0.5)) == 0
    end
end

function adaptive_attempt_allocs(problem, solver)
    cache=initialize_cache(problem, solver)
    sol=initialize_solution(problem, solver)
    function attempt!()
        cache.n=1
        cache.m=1
        NSDERungeKutta.step!(cache, sol, problem, solver, 0.01)
        NSDERungeKutta.adaptivestep!(cache, sol, solver, solver.adaptive, 0.01)
        return nothing
    end
    attempt!()
    return @allocated attempt!()
end
@testset "Adaptive attempt does not allocate solver settings" begin
    p=Logistic(0.3, (0.0, 1.0))
    @test adaptive_attempt_allocs(p, DormandPrince54(h=0.1)) == 0
end
