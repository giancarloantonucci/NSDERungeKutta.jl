struct ExplicitExponentialRungeKuttaSolver{stepsize_T} <: RungeKuttaSolver
    stepsize::stepsize_T
end

ExplicitExponentialRungeKuttaSolver(h::Real) = ExplicitExponentialRungeKuttaSolver(StepSize(h))
@doc (@doc ExplicitExponentialRungeKuttaSolver) EERK(args...; kwargs...) = ExplicitExponentialRungeKuttaSolver(args...; kwargs...)

mutable struct ExplicitExponentialRungeKuttaCache{n_T, Q_T, α_T, β_T, γ_T, E_T, E2_T} <: RungeKuttaCache
    n::n_T
    Q::Q_T
    α::α_T
    β::β_T
    γ::γ_T
    E::E_T
    E2::E2_T
end

function ExplicitExponentialRungeKuttaCache(problem::InitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver)
    @↓ L = problem.rhs
    @↓ h = solver.stepsize

    n = 1

    M = 16
    r = @. exp(π * 1im * ((1:M) - 0.5) / M)

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

function (solver::ExplicitExponentialRungeKuttaSolver)(problem::InitialValueProblem)
    solve(problem, solver)
end

function step!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver, cache::ExplicitExponentialRungeKuttaCache)
    @↓ n, Q, α, β, γ, E, E2 = cache
    @↓ u, t = solution
    @↓ f = problem.rhs.rhs
    @↓ h = solver.stepsize

    uₙ = u[n]
    tₙ = t[n+1] = t[n] + h

    Nu = f(uₙ, tₙ)
    aₙ = @. E2 * uₙ + Q * Nu
    Na = f(aₙ, tₙ)
    bₙ = @. E2 * uₙ + Q * Na
    Nb = f(bₙ, tₙ)
    cₙ = @. E2 * aₙ + Q * (2Nb - Nu)
    Nc = f(cₙ, tₙ)

    @. u[n+1] = E * uₙ + α * Nu + 2β * (Na + Nb) + γ * Nc
    return u[n+1], t[n+1]
end

adaptive_step!(solution::RungeKuttaSolution, solver::ExplicitExponentialRungeKuttaSolver, cache::RungeKuttaCache) = adaptive_step!(solution, solver, nothing, cache)

function ExponentialRK4(; h = 0.0)
    return EERK(h)
end
@doc (@doc ExponentialRK4) ERK4(args...; kwargs...) = ExponentialRK4(args...; kwargs...)
