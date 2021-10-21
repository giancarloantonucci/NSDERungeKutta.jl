"""
    step!(solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::ImplicitRungeKuttaSolver, cache::ImplicitRungeKuttaCache; savestages::Bool = false)

computes a step of the [`AbstractRungeKuttaSolution`](@ref) of an [`AbstractInitialValueProblem`](@ref) using an [`ImplicitRungeKuttaSolver`](@ref).
"""
function step!(cache::ImplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::ImplicitRungeKuttaSolver; savestages::Bool=false)
    @↓ n, v, Δk, J = cache
    @↓ u, t = solution
    k = savestages ? solution.k[n] : cache.k
    @↓ f!, Df! = problem.rhs
    @↓ s, A, c, b = solver.tableau
    @↓ ϵ, K = solver.newton
    @↓ h = solver.stepsize
    v = u[n+1] # avoid allocs
    # TO-DO: @← J = Df(v, u[n], t[n])
    Df!(J, v, u[n], t[n])
    Z = factorize(I - h * kron(A, J))
    # Compute stages
    zero!(k)
    for l = 1:K
        for i = 1:s
            zero!(v)
            for j in 1:s
                @. v += A[i, j] * k[j]
            end
            @. v = u[n] + h * v
            # TO-DO: @← Δk[i] = f(v, t[n] + h * c[i])
            f!(Δk[i], v, t[n] + h * c[i])
            @. Δk[i] -= k[i]
        end
        # TO-DO: ldiv!(Z, Δk)
        Δk_ = vcat(Δk...)
        ldiv!(Z, Δk_)
        if norm(Δk_) < ϵ * norm(k)
            break
        end
        # TO-DO: k .+= Δk
        for i in eachindex(k)
            L = length(k[i])
            k[i] .+= Δk_[((i - 1) * L + 1):(i * L)]
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
