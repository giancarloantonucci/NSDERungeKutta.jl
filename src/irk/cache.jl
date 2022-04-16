mutable struct ImplicitRungeKuttaCache{n_T<:Integer, m_T<:Integer, k_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), Δk_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), V_T<:(AbstractVector{ℂ} where ℂ<:Number), J_T<:(AbstractMatrix{ℂ} where ℂ<:Number), e_T<:(Ref{ℝ} where ℝ<:Real)} <: AbstractRungeKuttaCache
    n::n_T # step counter
    m::m_T # adaptive counter
    k::k_T # stages at step `n`
    Δk::Δk_T # stages update
    V::V_T # `vec(k)` or `vec(Δk)`
    J::J_T # Jacobian of RHS
    e::e_T # compensated summation error
end

function ImplicitRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitRungeKuttaSolver)
    n = m = 1
    @↓ u0_T ← eltype(u0), d ← length(u0) = problem
    @↓ s = solver.tableau
    k = Vector{u0_T}(undef, s, d)
    Δk = Vector{u0_T}(undef, s, d)
    V = Vector{u0_T}(undef, s*d)
    J = Matrix{u0_T}(undef, d, d)
    e = Ref(0.0)
    return ImplicitRungeKuttaCache(n, m, k, Δk, V, J, e)
end

#####
##### Functions
#####

RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitRungeKuttaSolver) = ImplicitRungeKuttaCache(problem, solver)
