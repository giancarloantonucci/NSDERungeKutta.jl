"""
    ExplicitRungeKuttaCache <: AbstractRungeKuttaCache

A composite type for the [`AbstractRungeKuttaCache`](@ref) of an [`ExplicitRungeKuttaSolver`](@ref).

# Constructors
```julia
ExplicitRungeKuttaCache(n, m, v, k)
ExplicitRungeKuttaCache(problem, solver)
```

# Arguments
- `n       :: Integer`                                               : step counter.
- `m       :: Integer`                                               : adaptive correction counter.
- `v       :: AbstractVector{Union{Number, AbstractVector{Number}}}` : temp for `solution.u[n]`.
- `k       :: AbstractVector{Union{Number, AbstractVector{Number}}}` : stages.
- `problem :: AbstractInitialValueProblem`                           : initial value problem..
- `solver  :: ExplicitRungeKuttaSolver`                              : explicit Runge-Kutta solver.
"""
mutable struct ExplicitRungeKuttaCache{n_T, m_T, v_T, k_T} <: AbstractRungeKuttaCache
    n::n_T
    m::m_T
    v::v_T
    k::k_T
end

function ExplicitRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u0)
    k = Vector{eltype(u0)}(undef, s, length(u0))
    return ExplicitRungeKuttaCache(n, m, v, k)
end

"""
    ImplicitRungeKuttaCache <: AbstractRungeKuttaCache

A composite type for the [`AbstractRungeKuttaCache`](@ref) of an [`ImplicitRungeKuttaSolver`](@ref).

# Constructors
```julia
ImplicitRungeKuttaCache(n, m, v, k, Δk, J)
ImplicitRungeKuttaCache(problem, solver)
```

# Arguments
- `n       :: Integer`                                               : step counter.
- `m       :: Integer`                                               : adaptive correction counter.
- `v       :: AbstractVector{Union{Number, AbstractVector{Number}}}` : temp for `solution.u[n]`.
- `k       :: AbstractVector{Union{Number, AbstractVector{Number}}}` : stages.
- `Δk      :: AbstractVector{Union{Number, AbstractVector{Number}}}` : stages' correction.
- `J       :: AbstractMatrix{Union{Number, AbstractMatrix{Number}}}` : Jacobian of RHS function.
- `problem :: AbstractInitialValueProblem`                           : initial value problem..
- `solver  :: ImplicitRungeKuttaSolver`                              : implicit Runge-Kutta solver.
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

mutable struct ExplicitExponentialRungeKuttaCache{n_T, Q_T, α_T, β_T, γ_T, E_T, E2_T} <: AbstractRungeKuttaCache
    n::n_T
    Q::Q_T
    α::α_T
    β::β_T
    γ::γ_T
    E::E_T
    E2::E2_T
end

function ExplicitExponentialRungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver)
    @↓ L = problem.rhs
    @↓ h = solver.stepsize
    #
    n = 1
    #
    M = 16
    r = @. exp(π * 1im * ((1:M) - 0.5) / M)
    #
    Z = h * L .+ r'
    Q = vec(h/M * real.(sum( @. (exp(Z/2) - 1.0) / Z                         ; dims = 2)))
    α = vec(h/M * real.(sum( @. (-4.0 - Z + exp(Z) * (4.0 - 3Z + Z^2)) / Z^3 ; dims = 2)))
    β = vec(h/M * real.(sum( @. (2.0 + Z + exp(Z) * (-2.0 + Z)) / Z^3        ; dims = 2)))
    γ = vec(h/M * real.(sum( @. (-4.0 - 3Z - Z^2 + exp(Z) * (4.0 - Z)) / Z^3 ; dims = 2)))
    # L ≡ diagonal matrix
    E = @. exp(h * L)
    E2 = @. exp(h/2 * L)
    return ExplicitExponentialRungeKuttaCache(n, Q, α, β, γ, E, E2)
end

RungeKuttaCache(problem, solver::ExplicitRungeKuttaSolver) = ExplicitRungeKuttaCache(problem, solver)
RungeKuttaCache(problem, solver::ImplicitRungeKuttaSolver) = ImplicitRungeKuttaCache(problem, solver)
RungeKuttaCache(problem, solver::ExplicitExponentialRungeKuttaSolver) = ExplicitExponentialRungeKuttaCache(problem, solver)
