mutable struct ImplicitRungeKuttaCache{n_T<:Integer, m_T<:Integer, v_T<:AbstractVector{<:Number}, k_T<:AbstractVector{<:AbstractVector{<:Number}}, Δk_T<:AbstractVector{<:AbstractVector{<:Number}}, V_T<:AbstractVector{<:Number}, J_T<:(AbstractMatrix{ℂ} where ℂ<:Number), e_T<:Ref{<:Real}} <: AbstractRungeKuttaCache
    n::n_T # step counter
    m::m_T # adaptive counter
    v::v_T # avoids allocation inside `adaptivestep!`
    k::k_T # stages at step `n`
    Δk::Δk_T # stages update
    V::V_T # `vec(k)` or `vec(Δk)`
    J::J_T # Jacobian of RHS
    e::e_T # compensated summation error
end

function ImplicitRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u0)
    k = [similar(u0) for i = 1:s]
    Δk = [similar(u0) for i = 1:s]
    d = length(u0)
    V = similar(u0, s*d)
    J = similar(u0, d, d)
    e = Ref(0.0)
    return ImplicitRungeKuttaCache(n, m, v, k, Δk, V, J, e)
end

#---------------------------------- FUNCTIONS ----------------------------------

RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitRungeKuttaSolver) = ImplicitRungeKuttaCache(problem, solver)
