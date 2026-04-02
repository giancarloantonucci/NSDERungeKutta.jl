# NSDERungeKutta/src/erk/cache.jl

mutable struct ExplicitRungeKuttaCache{
            n_T <: Integer,
            e_T <: Ref{<:AbstractFloat},
            v_T <: AbstractVector{<:Number},
            k_T <: AbstractVector{v_T}
        } <: AbstractRungeKuttaCache
    n :: n_T  # step counter
    m :: n_T  # adaptive step counter
    e :: e_T  # compensated summation error
    v :: v_T  # avoids allocation for `u[n+1]`
    k :: k_T  # stages at step `n`
end

function RungeKuttaCache(
            problem::AbstractInitialValueProblem,
            solver::ExplicitRungeKuttaSolver
        )
    @↓ u0, tspan = problem
    t0, tN = tspan
    @↓ s = solver.tableau
    n = m = 1
    e = Ref(zero(eltype(t0)))
    v = similar(u0)
    k = [similar(u0) for i = 1:s]
    return ExplicitRungeKuttaCache(n, m, e, v, k)
end
