"""
    RungeKuttaSolution(u, t) -> RungeKuttaSolution

returns a constructor for the `RungeKuttaSolution` of an `InitialValueProblem` with
- `u` : numerical solution.
- `t` : time grid.
"""
struct RungeKuttaSolution{u_T, t_T} <: InitialValueSolution
    u::u_T
    t::t_T
end

function RungeKuttaSolution(problem::InitialValueProblem, solver::RungeKuttaSolver)
    @↓ u0, tspan = problem
    @↓ h = solver
    T = eltype(u0)
    L = length(u0)
    t0, tN = tspan
    N = floor(Int, (tN - t0) / h) + 1
    u = BlockVector{T}(undef, [L for i = 1:N])
    u[Block(1)] = u0
    t = Vector{eltype(t0)}(undef, N)
    t[1] = t0
    return RungeKuttaSolution(u, t)
end

Base.length(solution::RungeKuttaSolution) = length(solution.t)
function Base.getindex(solution::RungeKuttaSolution, i::Int)
    @↓ u, t = solution
    N = length(t)
    return RungeKuttaSolution([u[BlockIndex(n, i)] for n = 1:N], t)
end
