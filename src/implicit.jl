"""
    ImplicitRungeKuttaSolver <: RungeKuttaSolver

A composite type for implicit [`RungeKuttaSolver`](@ref)s.

# Constructors
```julia
ImplicitRungeKuttaSolver(tableau, h, newton[, adaptive])
IRK(args...; kwargs...)
```

# Arguments
- `tableau :: ButcherTableau` : Butcher tableau.
- `stepsize :: StepSize` : step-size.
- `newton :: NewtonParameters` : simplified Newton's parameters.
- `adaptive :: AdaptiveParameters` : embedded method's parameters.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct ImplicitRungeKuttaSolver{tableau_T, stepsize_T, newton_T, adaptive_T} <: RungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    newton::newton_T
    adaptive::adaptive_T
end

ImplicitRungeKuttaSolver(tableau, h::Real, newton, adaptive) = ImplicitRungeKuttaSolver(tableau, StepSize(h), newton, adaptive)
ImplicitRungeKuttaSolver(tableau, stepsize, newton) = ImplicitRungeKuttaSolver(tableau, stepsize, newton, nothing)
@doc (@doc ImplicitRungeKuttaSolver) IRK(args...; kwargs...) = ImplicitRungeKuttaSolver(args...; kwargs...)

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

"""
    show(io::IO, solver::ImplicitRungeKuttaSolver)

prints a full description of `solver` and its contents to a stream `io`.
"""
Base.show(io::IO, solver::ImplicitRungeKuttaSolver) = NSDEBase._show(io, solver)

"""
    summary(io::IO, solver::ImplicitRungeKuttaSolver)

prints a brief description of `solver` to a stream `io`.
"""
Base.summary(io::IO, solver::ImplicitRungeKuttaSolver) = NSDEBase._summary(io, solver)

# ---------------------------------------------------------------------------- #
#                                    Methods                                   #
# ---------------------------------------------------------------------------- #

function (solver::ImplicitRungeKuttaSolver)(problem::InitialValueProblem; save_stages::Bool = false)
    solve(problem, solver; save_stages=save_stages)
end

"""
    ImplicitRungeKuttaCache <: RungeKuttaCache

A composite type for the [`RungeKuttaCache`](@ref) of an [`ImplicitRungeKuttaSolver`](@ref).

# Constructors
```julia
ImplicitRungeKuttaCache(n, m, v, k, Δk, J)
ImplicitRungeKuttaCache(problem, solver)
```

# Arguments
- `n  :: Integer` : step counter.
- `m  :: Integer` : adaptive correction counter.
- `v  :: AbstractVector{Union{Number, AbstractVector{Number}}}` : temp for `solution.u[n]`.
- `k  :: AbstractVector{Union{Number, AbstractVector{Number}}}` : stages.
- `Δk :: AbstractVector{Union{Number, AbstractVector{Number}}}` : stages' correction.
- `J  :: AbstractMatrix{Union{Number, AbstractMatrix{Number}}}` : Jacobian of RHS function.
- `problem :: InitialValueProblem`.
- `solver :: ImplicitRungeKuttaSolver`.
"""
mutable struct ImplicitRungeKuttaCache{n_T, m_T, v_T, k_T, Δk_T, J_T} <: RungeKuttaCache
    n::n_T
    m::m_T
    v::v_T
    k::k_T
    Δk::Δk_T
    J::J_T
end

function ImplicitRungeKuttaCache(problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u0)
    u0_T = eltype(u0)
    L = length(u0)
    k = Vector{u0_T}(undef, s, L)
    Δk = Vector{u0_T}(undef, s, L)
    J = Matrix{u0_T}(undef, L, L)
    return ImplicitRungeKuttaCache(n, m, v, k, Δk, J)
end

"""
    step!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver, cache::ImplicitRungeKuttaCache)

computes a step of the [`RungeKuttaSolution`](@ref) of an [`InitialValueProblem`](@ref) using an [`ImplicitRungeKuttaSolver`](@ref).
"""
function step!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver, cache::ImplicitRungeKuttaCache)
    @↓ n, v, Δk, J = cache
    @↓ u, t = solution
    k = solution.k isa Nothing ? cache.k : solution.k[n]
    @↓ f!, Df! = problem.rhs
    @↓ s, A, c, b = solver.tableau
    @↓ ϵ, K = solver.newton
    @↓ h = solver.stepsize
    v = u[n+1] # avoid allocs
    # @← J = Df(v, u[n], t[n])
    Df!(J, v, u[n], t[n])
    Z = factorize(I - h * kron(A, J))
    # Compute stages
    zero!(k)
    for l = 1:K
        for i = 1:s
            zero!(v)
            for j = 1:s
                @. v += A[i,j] * k[j]
            end
            @. v = u[n] + h * v
            # @← Δk[i] = f(v, t[n] + h * c[i])
            f!(Δk[i], v, t[n] + h * c[i])
            @. Δk[i] -= k[i]
        end
        # temporary workaround for ldiv!(Z, Δk)
        Δk_ = vcat(Δk...)
        ldiv!(Z, Δk_)
        if norm(Δk_) < ϵ * norm(k)
            break
        end
        # temporary workaround for k .+= Δk
        for i in eachindex(k)
            L = length(k[i])
            k[i] .+= Δk_[(i-1)*L+1:i*L]
        end
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
