mutable struct DiagonallyImplicitRungeKuttaCache{n_T<:Integer, e_T<:Ref{<:Real}, v_T<:AbstractVector{<:Number}, k_T<:AbstractVector{<:AbstractVector{<:Number}}, J_T<:AbstractMatrix{<:Number}} <: AbstractRungeKuttaCache
    n :: n_T # step counter
    m :: n_T # adaptive step counter
    e :: e_T # compensation error in `kahansum()`
    v :: v_T # avoids allocation for `u[n+1]`
    Uᵢ :: v_T # avoids allocating intermediate stages in `step!`
    Δkᵢ :: v_T # avoids allocating stages update in `step!`
    k :: k_T # stages at step `n`
    J :: J_T # avoids allocating Jacobian matrix of RHS in `step!`
end

function RungeKuttaCache(problem::AbstractInitialValueProblem, solver::DiagonallyImplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.tableau
    n = m = 1
    e = Ref(0.0)
    v = similar(u0)
    Uᵢ = similar(u0)
    Δkᵢ = similar(u0)
    k = [similar(u0) for i = 1:s]
    d = length(u0)
    J = similar(u0, d, d)
    return DiagonallyImplicitRungeKuttaCache(n, m, e, v, Uᵢ, Δkᵢ, k, J)
end
