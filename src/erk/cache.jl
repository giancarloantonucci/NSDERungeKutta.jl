mutable struct ExplicitRungeKuttaCache{n_T<:Integer, m_T<:Integer, v_T<:AbstractVector{<:Number}, k_T<:AbstractVector{<:AbstractVector{<:Number}}, e_T<:Ref{<:Real}} <: AbstractRungeKuttaCache
    n::n_T # step counter
    m::m_T # adaptive counter
    v::v_T # avoids allocation inside `adaptivestep!`
    k::k_T # stages at step `n`
    e::e_T # compensated summation error
end

function ExplicitRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u0)
    k = [similar(u0) for i = 1:s]
    e = Ref(0.0)
    return ExplicitRungeKuttaCache(n, m, v, k, e)
end

#---------------------------------- FUNCTIONS ----------------------------------

RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver) = ExplicitRungeKuttaCache(problem, solver)
