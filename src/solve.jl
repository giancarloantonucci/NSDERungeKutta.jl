"""
    solve!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver) :: AbstractRungeKuttaSolution

returns the [`AbstractRungeKuttaSolution`](@ref) of an [`AbstractInitialValueProblem`](@ref).
"""
function NSDEBase.solve!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver)
    cache = RungeKuttaCache(problem, solver)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ u, t = solution
    @↓ n = cache
    u0_T = eltype(u0)
    u0_L = length(u0)
    N = length(t)
    # WHILE instead of FOR loop -> adaptive methods
    while n < N && t[n] < tN
        step!(cache, solution, problem, solver)
        adaptivecheck!(cache, solution, solver)
        @↓ n = cache
        # Stop if update is too small -> adaptive methods
        if n > 1 && t[n] ≈ t[n - 1]
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
    solve(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver; savestages::Bool=false) :: AbstractRungeKuttaSolution

returns the [`AbstractRungeKuttaSolution`](@ref) of an [`AbstractInitialValueProblem`](@ref).
"""
function NSDEBase.solve(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver; savestages::Bool=false)
    solution = RungeKuttaSolution(problem, solver; savestages=savestages)
    solve!(solution, problem, solver)
    return solution
end
