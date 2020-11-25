mutable struct RungeKuttaCache{n_T, m_T, k_T, Δk_T, J_T}
    n::n_T
    m::m_T
    k::k_T
    Δk::Δk_T
    J::J_T
end

function RungeKuttaCache(problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver)
    @↓ u0, T ← eltype(u0), L ← length(u0) = problem
    @⤓ s = solver
    n = m = 1
    k = BlockVector{T}(undef, [L for i = 1:s])
    return RungeKuttaCache(n, m, k, ∅, ∅)
end

function RungeKuttaCache(problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver)
    @↓ u0, T ← eltype(u0), L ← length(u0) = problem
    @⤓ s = solver
    n = m = 1
    k = BlockVector{T}(undef, [L for i = 1:s])
    Δk = BlockVector{T}(undef, [L for i = 1:s])
    J = Matrix{T}(undef, L, L)
    return RungeKuttaCache(n, m, k, Δk, J)
end

function RungeKuttaCache(problem::InitialValueProblem, solver::AdditiveRungeKuttaSolver)
    @↓ u0, T ← eltype(u0), L ← length(u0) = problem
    @⤓ s = solver
    n = m = 1
    kᴵ = BlockVector{T}(undef, [L for i = 1:s])
    kᴱ = BlockVector{T}(undef, [L for i = 1:s])
    J = Matrix{T}(undef, L, L)
    return RungeKuttaCache(n, m, kᴵ, kᴱ, J)
end
