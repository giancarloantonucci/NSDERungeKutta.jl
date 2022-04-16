mutable struct DiagonallyImplicitRungeKuttaCache{n_T<:Integer, m_T<:Integer, k_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), Uᵢ_T<:(AbstractVector{ℂ} where ℂ<:Number), Δkᵢ_T<:(AbstractVector{ℂ} where ℂ<:Number), J_T<:(AbstractMatrix{ℂ} where ℂ<:Number), e_T<:(Ref{ℝ} where ℝ<:Real)} <: AbstractRungeKuttaCache
    n::n_T # step counter
    m::m_T # adaptive counter
    k::k_T # stages at step `n`
    Uᵢ::Uᵢ_T # ...
    Δkᵢ::Δkᵢ_T # stages update
    J::J_T # Jacobian of RHS
    e::e_T # compensated summation error
end

function DiagonallyImplicitRungeKuttaCache(problem::AbstractInitialValueProblem, solver::DiagonallyImplicitRungeKuttaSolver)
    n = m = 1
    @↓ u0 = problem
    u0_T = eltype(u0)
    @↓ s = solver.tableau
    d = length(u0)
    k = Vector{u0_T}(undef, s, d)
    Uᵢ = similar(u0)
    Δkᵢ = similar(u0)
    J = Matrix{u0_T}(undef, d, d)
    e = Ref(0.0)
    return DiagonallyImplicitRungeKuttaCache(n, m, k, Uᵢ, Δkᵢ, J, e)
end

#####
##### Functions
#####

RungeKuttaCache(problem::AbstractInitialValueProblem, solver::DiagonallyImplicitRungeKuttaSolver) = DiagonallyImplicitRungeKuttaCache(problem, solver)
