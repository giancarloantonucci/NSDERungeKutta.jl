mutable struct RungeKuttaCache{n_T, m_T, v_T, k_T, Δk_T, J_T}
    n::n_T
    m::m_T
    v::v_T
    k::k_T
    Δk::Δk_T
    J::J_T
end

function RungeKuttaCache(problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u0)
    T = eltype(u0)
    L = length(u0)
    k = BlockVector{T}(undef, [L for i = 1:s])
    return RungeKuttaCache(n, m, v, k, nothing, nothing)
end

function RungeKuttaCache(problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u0)
    T = eltype(u0)
    L = length(u0)
    k = BlockVector{T}(undef, [L for i = 1:s])
    Δk = BlockVector{T}(undef, [L for i = 1:s])
    J = Matrix{T}(undef, L, L)
    return RungeKuttaCache(n, m, v, k, Δk, J)
end
