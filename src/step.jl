"""
step!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver, cache::RungeKuttaCache)

computes a step of the `RungeKuttaSolution` of an `InitialValueProblem` using an `ExplicitRungeKuttaSolver`.
"""
function step!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver, cache::RungeKuttaCache, save_stages)
    @↓ n, v = cache
    @↓ u, t = solution
    k = if save_stages
        @↓ k ← k[n] = solution
    else
        @↓ k = cache
    end
    @↓ f! = problem.rhs
    @↓ tableau, h = solver
    @↓ s, A, c, b = tableau
    v = u[n+1] # avoid allocs
    # compute stages
    for i = 1:s
        zero!(v)
        for j = 1:i-1
            @. v += A[i,j] * k[j]
        end
        @. v = u[n] + h * v
        # @← k[i] = f(v, t[n] + h * c[i])
        f!(k[i], v, t[n] + h * c[i])
        t[n+1] = t[n] + h
    end
    # compute step
    zero!(v)
    for i = 1:s
        @. v += b[i] * k[i]
    end
    @. v = u[n] + h * v
    t[n+1] = t[n] + h
end

"""
step!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver, cache::RungeKuttaCache)

computes a step of the `RungeKuttaSolution` of an `InitialValueProblem` using an `ImplicitRungeKuttaSolver`.
"""
function step!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver, cache::RungeKuttaCache, save_stages)
    @↓ n, v, Δk, J = cache
    @↓ u, t = solution
    k = if save_stages
        @↓ k ← k[n] = solution
    else
        @↓ k = cache
    end
    @↓ f!, Df! = problem.rhs
    @↓ tableau, h, ϵ, K = solver
    @↓ s, A, c, b = tableau
    v = u[n+1] # avoid allocs
    # @← J = Df(v, u[n], t[n])
    Df!(J, v, u[n], t[n])
    Z = factorize(I - h * kron(A, J))
    # compute stages
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
    # compute step
    zero!(v)
    for i = 1:s
        @. v += b[i] * k[i]
    end
    @. v = u[n] + h * v
    t[n+1] = t[n] + h
end
