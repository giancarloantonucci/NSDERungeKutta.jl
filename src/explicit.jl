"""
    ExplicitRungeKuttaSolver <: RungeKuttaSolver

A composite type for explicit [`RungeKuttaSolver`](@ref)s.

# Constructors
```julia
ExplicitRungeKuttaSolver(tableau, h[, adaptive])
ERK(args...; kwargs...)
```

# Arguments
- `tableau :: ButcherTableau` : Butcher tableau.
- `stepsize :: StepSize` : step-size.
- `adaptive :: AdaptiveParameters` : embedded method's parameters.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct ExplicitRungeKuttaSolver{tableau_T, stepsize_T, adaptive_T} <: RungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    adaptive::adaptive_T
end

ExplicitRungeKuttaSolver(tableau, h::Real, adaptive) = ExplicitRungeKuttaSolver(tableau, StepSize(h), adaptive)
ExplicitRungeKuttaSolver(tableau, stepsize) = ExplicitRungeKuttaSolver(tableau, stepsize, nothing)
@doc (@doc ExplicitRungeKuttaSolver) ERK(args...; kwargs...) = ExplicitRungeKuttaSolver(args...; kwargs...)

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

"""
    show(io::IO, solver::ExplicitRungeKuttaSolver)

prints a full description of `solver` and its contents to a stream `io`.
"""
Base.show(io::IO, solver::ExplicitRungeKuttaSolver) = _show(io, solver)

"""
    summary(io::IO, solver::ExplicitRungeKuttaSolver)

prints a brief description of `solver` to a stream `io`.
"""
Base.summary(io::IO, solver::ExplicitRungeKuttaSolver) = _summary(io, solver)

# ---------------------------------------------------------------------------- #
#                                    Methods                                   #
# ---------------------------------------------------------------------------- #

function (solver::ExplicitRungeKuttaSolver)(problem::InitialValueProblem; save_stages::Bool = false)
    solve(problem, solver; save_stages=save_stages)
end

"""
    ExplicitRungeKuttaCache <: RungeKuttaCache

A composite type for the [`RungeKuttaCache`](@ref) of an [`ExplicitRungeKuttaSolver`](@ref).

# Constructors
```julia
ExplicitRungeKuttaCache(n, m, v, k)
ExplicitRungeKuttaCache(problem, solver)
```

# Arguments
- `n  :: Integer` : step counter.
- `m  :: Integer` : adaptive correction counter.
- `v  :: AbstractVector{Union{Number, AbstractVector{Number}}}` : temp for `solution.u[n]`.
- `k  :: AbstractVector{Union{Number, AbstractVector{Number}}}` : stages.
- `problem :: InitialValueProblem`.
- `solver :: ExplicitRungeKuttaSolver`.
"""
mutable struct ExplicitRungeKuttaCache{n_T, m_T, v_T, k_T} <: RungeKuttaCache
    n::n_T
    m::m_T
    v::v_T
    k::k_T
end

function ExplicitRungeKuttaCache(problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u0)
    k = Vector{eltype(u0)}(undef, s, length(u0))
    return ExplicitRungeKuttaCache(n, m, v, k)
end

"""
    step!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver, cache::ExplicitRungeKuttaCache)

computes a step of the [`RungeKuttaSolution`](@ref) of an [`InitialValueProblem`](@ref) using an [`ExplicitRungeKuttaSolver`](@ref).
"""
function step!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver, cache::ExplicitRungeKuttaCache)
    @↓ n, v = cache
    @↓ u, t = solution
    k = solution.k isa Nothing ? cache.k : solution.k[n]
    @↓ f! = problem.rhs
    @↓ s, A, c, b = solver.tableau
    @↓ h = solver.stepsize
    v = u[n+1] # avoid allocs
    # Compute stages
    for i = 1:s
        zero!(v)
        for j = 1:i-1
            @. v += A[i,j] * k[j]
        end
        @. v = u[n] + h * v
        # @← k[i] = f(v, t[n] + h * c[i])
        f!(k[i], v, t[n] + h * c[i])
    end
    # Compute step
    zero!(v)
    for i = 1:s
        @. v += b[i] * k[i]
    end
    @. v = u[n] + h * v
    t[n+1] = t[n] + h
    return u[n+1], t[n+1]
end
