
"""
    ImplicitRungeKuttaCache <: AbstractRungeKuttaCache

A composite type for the [`AbstractRungeKuttaCache`](@ref) of an [`ImplicitRungeKuttaSolver`](@ref).

# Constructors
```julia
ImplicitRungeKuttaCache(n, m, v, k, Δk, J)
ImplicitRungeKuttaCache(problem, solver)
```

## Arguments
- `n :: Integer` : step counter.
- `m :: Integer` : adaptive correction counter.
- `v :: AbstractVector{<:Number}` : temp for `solution.u[n]`.
- `k :: AbstractVector{<:Union{Number, AbstractVector{<:Number}}}` : stages.
- `Δk :: AbstractVector{<:Union{Number, AbstractVector{<:Number}}}` : stages' correction.
- `J :: AbstractMatrix{<:Number}` : Jacobian of right-hand side derivative.
- `problem :: AbstractInitialValueProblem`
- `solver :: ExplicitRungeKuttaSolver`

# Functions
- [`RungeKuttaCache`](@ref) : alternative caller.
"""
mutable struct ImplicitRungeKuttaCache{n_T, m_T, v_T, k_T, Δk_T, J_T} <: AbstractRungeKuttaCache
    n::n_T
    m::m_T
    v::v_T
    k::k_T
    Δk::Δk_T
    J::J_T
end

function ImplicitRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitRungeKuttaSolver)
    n = m = 1
    @↓ u0 = problem
    v = similar(u0)
    @↓ s = solver.tableau
    u0_T = eltype(u0)
    u0_L = length(u0)
    k = Vector{u0_T}(undef, s, u0_L)
    Δk = Vector{u0_T}(undef, s, u0_L)
    J = Matrix{u0_T}(undef, u0_L, u0_L)
    return ImplicitRungeKuttaCache(n, m, v, k, Δk, J)
end

#####
##### Functions
#####

RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitRungeKuttaSolver) = ImplicitRungeKuttaCache(problem, solver)
