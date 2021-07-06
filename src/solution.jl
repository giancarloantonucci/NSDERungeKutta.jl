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
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ h = solver.stepsize
    @↓ s = solver.tableau
    N = floor(Int, (tN - t0) / h) + 1
    u = Vector{eltype(u0)}(undef, N, length(u0))
    u[1] = u0
    t = Vector{typeof(t0)}(undef, N)
    t[1] = t0
    k = save_stages ? Vector{uT}(undef, N, s, L) : nothing
    return RungeKuttaSolution(u, t, k)
end

Base.length(solution::RungeKuttaSolution) = length(solution.t)
