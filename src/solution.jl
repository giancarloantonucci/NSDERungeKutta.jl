"""
    RungeKuttaSolution(u, t) -> RungeKuttaSolution

returns a constructor for the numerical solution of an `InitialValueProblem`.

# Arguments
- `u` : numerical solution.
- `t` : time grid.
"""
struct RungeKuttaSolution{u_T, t_T} <: InitialValueSolution
    u::u_T
    t::t_T
end

function RungeKuttaSolution(problem::InitialValueProblem, solver::RungeKuttaSolver)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ h = solver
    # build t vector
    N = floor(Int, (tN - t0) / h) + 1
    tT = eltype(t0)
    t = Vector{tT}(undef, N)
    t[1] = t0
    # build u vector
    L = length(u0)
    uT = eltype(u0)
    u = BlockVector{uT}(undef, [L for i = 1:N])
    u[Block(1)] = u0
    return RungeKuttaSolution(u, t)
end

Base.length(solution::RungeKuttaSolution) = length(solution.t)
function Base.getindex(solution::RungeKuttaSolution, i::Int)
    @↓ u, t = solution
    N = length(t)
    u_i = [u[BlockIndex(n, i)] for n = 1:N] # change to BlockVector
    return RungeKuttaSolution(u_i, t)
end
