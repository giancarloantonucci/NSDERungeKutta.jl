mutable struct ImplicitExplicitRungeKuttaCache{n_T<:Integer, m_T<:Integer, kᴵ_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), kᴱ_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), Uᵢ_T<:(AbstractVector{ℂ} where ℂ<:Number), J_T<:(AbstractMatrix{ℂ} where ℂ<:Number), e_T<:(Ref{ℝ} where ℝ<:Real)} <: AbstractRungeKuttaCache
    n::n_T # step counter.
    m::m_T # adaptive counter.
    kᴵ::kᴵ_T # implicit stages at step `n`
    kᴱ::kᴱ_T # explicit stages at step `n`
    Uᵢ::Uᵢ_T # ...
    J::J_T # Jacobian of stiff part of RHS
    e::e_T # compensated summation error
end

function ImplicitExplicitRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitExplicitRungeKuttaSolver)
    n = m = 1
    @↓ u0_T ← eltype(u0), d ← length(u0) = problem
    @↓ s = solver.implicitableau
    kᴵ = Vector{u0_T}(undef, s, d)
    kᴱ = Vector{u0_T}(undef, s, d)
    Uᵢ = Vector{u0_T}(undef, d)
    J = Matrix{u0_T}(undef, d, d)
    e = Ref(0.0)
    return ImplicitExplicitRungeKuttaCache(n, m, kᴵ, kᴱ, Uᵢ, J, e)
end

#####
##### Functions
#####

RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitExplicitRungeKuttaSolver) = ImplicitExplicitRungeKuttaCache(problem, solver)
