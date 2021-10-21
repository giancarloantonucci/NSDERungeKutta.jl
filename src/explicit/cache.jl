"""
    ExplicitRungeKuttaCache <: AbstractRungeKuttaCache

A composite type for the [`AbstractRungeKuttaCache`](@ref) of an [`ExplicitRungeKuttaSolver`](@ref).

# Constructors
```julia
ExplicitRungeKuttaCache(n, m, v, k)
ExplicitRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver)
```

# Arguments
- `n :: Integer` : step counter.
- `m :: Integer` : adaptive correction counter.
- `v :: AbstractVector{Union{Number, AbstractVector{Number}}}` : temp for `solution.u[n]`.
- `k :: AbstractVector{Union{Number, AbstractVector{Number}}}` : stages.
"""
mutable struct ExplicitRungeKuttaCache{n_T, m_T, v_T, k_T} <: AbstractRungeKuttaCache
    n::n_T
    m::m_T
    v::v_T
    k::k_T
end

function ExplicitRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver)
    n = m = 1
    @↓ u0 = problem
    v = similar(u0)
    @↓ s = solver.tableau
    u0_T = eltype(u0)
    u0_L = length(u0)
    k = Vector{u0_T}(undef, s, u0_L)
    return ExplicitRungeKuttaCache(n, m, v, k)
end

#####
##### Functions
#####

function RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver)
    return ExplicitRungeKuttaCache(problem, solver)
end
