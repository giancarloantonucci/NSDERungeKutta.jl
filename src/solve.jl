Cache(problem, solver::ExplicitRungeKuttaSolver) = ExplicitRungeKuttaCache(problem, solver)
Cache(problem, solver::ImplicitRungeKuttaSolver) = ImplicitRungeKuttaCache(problem, solver)

"""
    solve!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::RungeKuttaSolver)

returns the `RungeKuttaSolution` of an `InitialValueProblem`.
"""
function NSDEBase.solve!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::RungeKuttaSolver)
    cache = Cache(problem, solver)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ u, t = solution
    @↓ n = cache
    u0_T = eltype(u0)
    u0_L = length(u0)
    # WHILE instead of FOR loop -> adaptive methods
    N = length(t)
    while n < N && t[n] < tN
        step!(solution, problem, solver, cache)
        adaptive_step!(solution, solver, cache)
        @↓ n = cache
        # Stop if update is too small -> adaptive methods
        if n > 1 && t[n] ≈ t[n-1]
            break
        end
        # Increase size of output vectors -> adaptive methods
        # if n == N && t[n] < tN
        #     N += 1
        #     push!(u, Vector{u0_T}(undef, u0_L))
        #     resize!(t, N)
        # end
    end
    # Decrease size of output vectors -> adaptive methods
    if n < N
        resize!(u, n)
        resize!(t, n)
    end
    return solution
end

"""
    solve(problem::InitialValueProblem, solver::RungeKuttaSolver[; save_stages::Bool = false) :: RungeKuttaSolution

returns the `RungeKuttaSolution` of an `InitialValueProblem`.
"""
function NSDEBase.solve(problem::InitialValueProblem, solver::RungeKuttaSolver; save_stages::Bool = false)
    solution = RungeKuttaSolution(problem, solver, save_stages)
    NSDEBase.solve!(solution, problem, solver)
    return solution
end
