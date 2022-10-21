mutable struct ExplicitExponentialRungeKuttaCache{n_T<:Integer, m_T<:Integer, v_T<:(AbstractVector{ℂ} where ℂ<:Number), k_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), e_T<:(Ref{ℝ} where ℝ<:Real)} <: AbstractRungeKuttaCache
    n::n_T # step counter
    m::m_T # adaptive counter
    v::v_T # avoids allocation for `u[n+1]`
    k::k_T # stages at step `n`
    e::e_T # compensated summation error
end

function ExplicitExponentialRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u0)
    k = [similar(u0) for i = 1:s]
    e = Ref(0.0)
    return ExplicitExponentialRungeKuttaCache(n, m, v, k, e)
end

#---------------------------------- FUNCTIONS ----------------------------------

RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver) = ExplicitExponentialRungeKuttaCache(problem, solver)
