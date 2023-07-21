mutable struct ImplicitExplicitRungeKuttaCache{n_T<:Integer, m_T<:Integer, v_T<:AbstractVector{<:Number}, kᴵ_T<:AbstractVector{<:AbstractVector{<:Number}}, kᴱ_T<:AbstractVector{<:AbstractVector{<:Number}}, Uᵢ_T<:AbstractVector{<:Number}, J_T<:(AbstractMatrix{ℂ} where ℂ<:Number), e_T<:Ref{<:Real}} <: AbstractRungeKuttaCache
    n::n_T # step counter
    m::m_T # adaptive counter
    v::v_T # avoids allocation inside `adaptivestep!`
    kᴵ::kᴵ_T # implicit stages at step `n`
    kᴱ::kᴱ_T # explicit stages at step `n`
    Uᵢ::Uᵢ_T # avoids allocation inside `step!`
    J::J_T # Jacobian of stiff part of RHS
    e::e_T # compensated summation error
end

function ImplicitExplicitRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitExplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.implicitableau
    n = m = 1
    kᴵ = [similar(u0) for i = 1:s]
    kᴱ = [similar(u0) for i = 1:s]
    d = length(u0)
    Uᵢ = similar(u0, d)
    J = similar(u0, d, d)
    e = Ref(0.0)
    return ImplicitExplicitRungeKuttaCache(n, m, v, kᴵ, kᴱ, Uᵢ, J, e)
end

#---------------------------------- FUNCTIONS ----------------------------------

RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitExplicitRungeKuttaSolver) = ImplicitExplicitRungeKuttaCache(problem, solver)
