mutable struct DiagonallyImplicitRungeKuttaCache{n_T<:Integer, m_T<:Integer, v_T<:(AbstractVector{ℂ} where ℂ<:Number), k_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), Uᵢ_T<:(AbstractVector{ℂ} where ℂ<:Number), Δkᵢ_T<:(AbstractVector{ℂ} where ℂ<:Number), J_T<:(AbstractMatrix{ℂ} where ℂ<:Number), e_T<:(Ref{ℝ} where ℝ<:Real)} <: AbstractRungeKuttaCache
    n::n_T # step counter
    m::m_T # adaptive counter
    v::v_T # avoids allocation inside `adaptivestep!`
    k::k_T # stages at step `n`
    Uᵢ::Uᵢ_T # avoids allocation inside `step!`
    Δkᵢ::Δkᵢ_T # avoids allocating `Δk`
    J::J_T # Jacobian of RHS
    e::e_T # compensated summation error
end

function DiagonallyImplicitRungeKuttaCache(problem::AbstractInitialValueProblem, solver::DiagonallyImplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u0)
    k = [similar(u0) for i = 1:s]
    Uᵢ = similar(u0)
    Δkᵢ = similar(u0)
    d = length(u0)
    J = similar(u0, d, d)
    e = Ref(0.0)
    return DiagonallyImplicitRungeKuttaCache(n, m, v, k, Uᵢ, Δkᵢ, J, e)
end

#---------------------------------- FUNCTIONS ----------------------------------

RungeKuttaCache(problem::AbstractInitialValueProblem, solver::DiagonallyImplicitRungeKuttaSolver) = DiagonallyImplicitRungeKuttaCache(problem, solver)
