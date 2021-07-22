"""
    Cache(n, m, v, k, Δk, J) <: RungeKuttaParameters

returns a constructor containing the temporary objects needed by a `RungeKuttaSolver`.

# Arguments
- `n  :: Integer`                                               : step counter.
- `m  :: Integer`                                               : adaptive correction counter.
- `v  :: AbstractVector{Union{Number, AbstractVector{Number}}}` : `solution.u[n]`.
- `k  :: AbstractVector{Union{Number, AbstractVector{Number}}}` : stages.
- `Δk :: AbstractVector{Union{Number, AbstractVector{Number}}}` : stages' correction.
- `J  :: AbstractMatrix{Union{Number, AbstractMatrix{Number}}}` : Jacobian of RHS function.
"""
mutable struct Cache{n_T, m_T, v_T, k_T, Δk_T, J_T}
    n::n_T
    m::m_T
    v::v_T
    k::k_T
    Δk::Δk_T
    J::J_T
end

Cache(n, m, v, k) = Cache(n, m, v, k, nothing, nothing)

function Cache(problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver)
    @↓ u₀ = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u₀)
    k = Vector{eltype(u₀)}(undef, s, length(u₀))
    return Cache(n, m, v, k)
end

function Cache(problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver)
    @↓ u₀ = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u₀)
    u₀_T = eltype(u₀)
    L = length(u₀)
    k = Vector{u₀_T}(undef, s, L)
    Δk = Vector{u₀_T}(undef, s, L)
    J = Matrix{u₀_T}(undef, L, L)
    return Cache(n, m, v, k, Δk, J)
end
