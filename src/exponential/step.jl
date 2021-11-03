"""
    step!(cache::ExplicitExponentialRungeKuttaCache, solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver)

computes a step of the [`AbstractRungeKuttaSolution`](@ref) of an [`AbstractInitialValueProblem`](@ref) using an [`ExplicitExponentialRungeKuttaSolver`](@ref).
"""
function step!(cache::ExplicitExponentialRungeKuttaCache, solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver; kwargs...)
    @↓ n, Q, α, β, γ, E, E2 = cache
    @↓ u, t = solution
    @↓ f = problem.rhs.rhs
    @↓ h = solver.stepsize

    uₙ = u[n]
    tₙ = t[n + 1] = t[n] + h

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
