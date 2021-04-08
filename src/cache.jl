mutable struct RungeKuttaCache{n_T, m_T, v_T, k_T, Δk_T, J_T}
    n::n_T
    m::m_T
    v::v_T
    k::k_T
    Δk::Δk_T
    J::J_T
end

RungeKuttaCache(n, m, v, k) = RungeKuttaCache(n, m, v, k, nothing, nothing)

function RungeKuttaCache(problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver)
    @↓ u0, T ← eltype(u0), L ← length(u0) = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u0)
    k = Vector{T}(undef, s, L)
    return RungeKuttaCache(n, m, v, k)
end

function RungeKuttaCache(problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver)
    @↓ u0, T ← eltype(u0), L ← length(u0) = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u0)
    k = Vector{T}(undef, s, L)
    Δk = Vector{T}(undef, s, L)
    J = Matrix{T}(undef, L, L)
    return RungeKuttaCache(n, m, v, k, Δk, J)
end
