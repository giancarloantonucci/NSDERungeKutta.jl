mutable struct ImplicitExplicitRungeKuttaCache{n_T <: Integer, e_T <: Ref{<:Real}, v_T <: AbstractVector{<:Number}, k_T <: AbstractVector{<:AbstractVector{<:Number}}, J_T <: AbstractMatrix{<:Number}, } <: AbstractRungeKuttaCache
    n :: n_T # step counter
    m :: n_T # adaptive step counter
    e :: e_T # compensated summation error
    v :: v_T # avoids allocation for `u[n+1]`
    Uᵢ :: v_T # avoids allocation inside `step!`
    kᴵ :: k_T # implicit stages at step `n`
    kᴱ :: k_T # explicit stages at step `n`
    J :: J_T # Jacobian of stiff part of RHS
end

function RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitExplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.implicitableau
    n = m = 1
    e = Ref(0.0)
    v = similar(u0)
    Uᵢ = similar(u0)
    kᴵ = [similar(u0) for i = 1:s]
    kᴱ = [similar(u0) for i = 1:s]
    d = length(u0)
    J = similar(u0, d, d)
    return ImplicitExplicitRungeKuttaCache(n, m, e, v, Uᵢ, kᴵ, kᴱ, J)
end
