"""
    step!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver, cache::ExplicitRungeKuttaCache)

computes a step of the [`AbstractRungeKuttaSolution`](@ref) of an [`AbstractInitialValueProblem`](@ref) using an [`ExplicitRungeKuttaSolver`](@ref).
"""
function step!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver, cache::ExplicitRungeKuttaCache)
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

"""
    step!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::ImplicitRungeKuttaSolver, cache::ImplicitRungeKuttaCache)

computes a step of the [`AbstractRungeKuttaSolution`](@ref) of an [`AbstractInitialValueProblem`](@ref) using an [`ImplicitRungeKuttaSolver`](@ref).
"""
function step!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::ImplicitRungeKuttaSolver, cache::ImplicitRungeKuttaCache)
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

function step!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::ExplicitExponentialRungeKuttaSolver, cache::ExplicitExponentialRungeKuttaCache)
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
