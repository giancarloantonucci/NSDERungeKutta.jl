mutable struct ExplicitExponentialRungeKuttaCache{n_T<:Integer, m_T<:Integer, v_T<:(AbstractVector{ℂ} where ℂ<:Number), k_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), e_T<:(Ref{ℝ} where ℝ<:Real)} <: AbstractRungeKuttaCache
    n::n_T
    m::m_T
    v::v_T
    k::k_T
    e::e_T
end

function ExplicitExponentialRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver)
    n = m = 1
    @↓ u0 = problem
    v = similar(u0)
    u0_T = eltype(u0)
    @↓ s = solver.tableau
    d = length(u0)
    k = Vector{u0_T}(undef, s, d)
    e = Ref(0.0)
    return ExplicitExponentialRungeKuttaCache(n, m, v, k, e)
end

#####
##### Functions
#####

RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver) = ExplicitExponentialRungeKuttaCache(problem, solver)
