"""
    ExplicitRungeKuttaCache <: AbstractRungeKuttaCache

A composite type for the [`AbstractRungeKuttaCache`](@ref) of an [`ExplicitRungeKuttaSolver`](@ref).

# Constructors
```julia
ExplicitRungeKuttaCache(n, m, v, k)
ExplicitRungeKuttaCache(problem, solver)
```

## Arguments
- `n :: Integer` : step counter.
- `m :: Integer` : adaptive correction counter.
- `v :: AbstractVector{<:Number}` : temp for `solution.u[n]`.
- `k :: AbstractVector{<:AbstractVector{<:Number}}` : stages.
- `problem :: AbstractInitialValueProblem`
- `solver :: ExplicitRungeKuttaSolver`

# Functions
- [`RungeKuttaCache`](@ref) : alternative caller.
"""
mutable struct ExplicitRungeKuttaCache{n_T<:Integer, m_T<:Integer, v_T<:AbstractVector{<:Number}, k_T<:AbstractVector{<:AbstractVector{<:Number}}} <: AbstractRungeKuttaCache
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
    k = Vector{eltype(u0)}(undef, s, length(u0))
    return ExplicitRungeKuttaCache(n, m, v, k)
end

#####
##### Functions
#####

"""
    RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver) :: ExplicitRungeKuttaCache
    RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver) :: ExplicitExponentialRungeKuttaCache
    RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ImplicitRungeKuttaSolver) :: ImplicitRungeKuttaCache
    
returns an [`AbstractRungeKuttaCache`](@ref) for each type of `solver`.
"""
RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver) = ExplicitRungeKuttaCache(problem, solver)
