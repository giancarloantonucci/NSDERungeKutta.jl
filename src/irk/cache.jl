mutable struct ImplicitRungeKuttaCache{n_T<:Integer, e_T<:Ref{<:Real}, v_T<:AbstractVector{<:Number}, k_T<:AbstractVector{<:AbstractVector{<:Number}}, J_T<:(AbstractMatrix{ℂ} where ℂ<:Number)} <: AbstractRungeKuttaCache
    n::n_T # step counter
    m::n_T # adaptive step counter
    e::e_T # compensated summation error
    v::v_T # avoids allocation for `u[n+1]`
    V::v_T # avoids allocation for `vec(k)` and `vec(Δk)`
    k::k_T # stages at step `n`
    Δk::k_T # stages update
    J::J_T # Jacobian of RHS
end

function RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.tableau
    n = m = 1
    e = Ref(0.0)
    v = similar(u0)
    d = length(u0)
    V = similar(u0, s*d)
    k = [similar(u0) for i = 1:s]
    Δk = [similar(u0) for i = 1:s]
    J = similar(u0, d, d)
    return ImplicitRungeKuttaCache(n, m, e, v, V, k, Δk, J)
end
