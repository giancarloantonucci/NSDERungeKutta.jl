"""
    RungeKuttaSolution(u, t[, k]) -> RungeKuttaSolution

returns a constructor for the numerical solution of an `InitialValueProblem`.

# Arguments
- `u` : numerical solution.
- `t` : time grid.
- `k` : stages.
"""
struct RungeKuttaSolution{u_T, t_T, k_T} <: InitialValueSolution
    u::u_T
    t::t_T
    k::k_T
end

RungeKuttaSolution(u, t) = RungeKuttaSolution(u, t, nothing)

function RungeKuttaSolution(problem::InitialValueProblem, solver::RungeKuttaSolver, save_stages)
    @↓ u0, uT ← eltype(u0), L ← length(u0) = problem
    @↓ (t0, tN) ← tspan, tT ← eltype(tspan) = problem
    @↓ h = solver
    @↓ s = solver.tableau
    N = floor(Int, (tN - t0) / h) + 1
    u = Vector{uT}(undef, N, L); u[1] = u0
    t = Vector{tT}(undef, N); t[1] = t0
    k = if save_stages
        Vector{uT}(undef, N, s, L)
    else
        nothing
    end
    return RungeKuttaSolution(u, t, k)
end

Base.length(solution::RungeKuttaSolution) = length(solution.t)
