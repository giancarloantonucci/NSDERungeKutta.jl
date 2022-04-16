mutable struct ExplicitRungeKuttaCache{n_T<:Integer, m_T<:Integer, v_T<:(AbstractVector{ℂ} where ℂ<:Number), k_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), e_T<:(Ref{ℝ} where ℝ<:Real)} <: AbstractRungeKuttaCache
    n::n_T # step counter
    m::m_T # adaptive counter
    v::v_T # ...
    k::k_T # stages at step `n`
    e::e_T # compensated summation error
end

function ExplicitRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver)
    n = m = 1
    @↓ u0 = problem
    v = similar(u0)
    u0_T = eltype(u0)
    @↓ s = solver.tableau
    d = length(u0)
    k = Vector{u0_T}(undef, s, d)
    e = Ref(0.0)
    return ExplicitRungeKuttaCache(n, m, v, k, e)
end

#####
##### Functions
#####

RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver) = ExplicitRungeKuttaCache(problem, solver)
