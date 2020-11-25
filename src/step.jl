"""
    step!(cache::RungeKuttaCache, solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver)

computes a step of the `RungeKuttaSolution` of an `InitialValueProblem` using an `ExplicitRungeKuttaSolver`.
"""
function step!(cache::RungeKuttaCache, solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver)
    @↓ n, k = cache
    @↓ u, t = solution
    @↓ f = problem.rhs
    @↓ tableau, h = solver
    @↓ A, b, c, s = tableau
    v = u[Block(n+1)] # to avoid allocs
    # compute stages
    for i = 1:s
        zero!(v)
        for j = 1:i-1
            @. v += A[i,j] * k[Block(j)]
        end
        @. v = u[Block(n)] + h * v
        @← k[Block(i)] = f(v, t[n] + h * c[i])
    end
    # compute step
    zero!(v)
    for i = 1:s
        @. v += b[i] * k[Block(i)]
    end
    @. v = u[Block(n)] + h * v
    t[n+1] = t[n] + h
end

"""
    step!(cache::RungeKuttaCache, solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver)

computes a step of the `RungeKuttaSolution` of an `InitialValueProblem` using an `ImplicitRungeKuttaSolver`.
"""
function step!(cache::RungeKuttaCache, solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver)
    @↓ n, k, Δk, J = cache
    @↓ u, t = solution
    @↓ f, Df = problem.rhs
    @↓ tableau, h, ϵ, K = solver
    @↓ A, b, c, s = tableau
    v = u[Block(n+1)] # to avoid allocs
    @← J = Df(u[Block(n)], t[n])
    Z = factorize(I - h * kron(A, J))
    zero!(k)
    # compute stages
    for l = 1:K
        for i = 1:s
            zero!(v)
            for j = 1:s
                @. v += A[i,j] * k[Block(j)]
            end
            @. v = u[Block(n)] + h * v
            @← Δk[Block(i)] = f(v, t[n] + h * c[i])
            @. Δk[Block(i)] -= k[Block(i)]
        end
        Δk .= Z \ Δk
        if norm(Δk) < ϵ * norm(k)
            break
        end
        k .+= Δk
    end
    # compute step
    zero!(v)
    for i = 1:s
        @. v += b[i] * k[Block(i)]
    end
    @. v = u[Block(n)] + h * v
    t[n+1] = t[n] + h
end

"""
    step!(cache::RungeKuttaCache, solution::RungeKuttaSolution, problem::InitialValueProblem, solver::AdditiveRungeKuttaSolver)

computes a step of the `RungeKuttaSolution` of an `InitialValueProblem` using an `AdditiveRungeKuttaSolver`.
"""
function step!(cache::RungeKuttaCache, solution::RungeKuttaSolution, problem::InitialValueProblem, solver::AdditiveRungeKuttaSolver)
    @↓ cache    = kᴱ, kᴵ, J
    @↓ solution = u, t
    @↓ problem  = f, g
    @↓ solver   = h, ϵᵣ, K, Aᴱ, bᴱ, cᴱ, Aᴵ, bᴵ, cᴵ, s
    v = u[n+1] # to avoid allocs
    ṽ   = similar(v)
    Δṽ  = similar(v)
    f_ṽ = similar(v)
    @← J = Df(u[n], t[n])
    for i = 1:s
        zero!(v)
        for j = 1:i-1
            @. v += Aᴵ[i,j] * kᴵ[j] + Aᴱ[i,j] * kᴱ[j]
        end
        @. v = u[n] + h * v
        Z = factorize(I - h * Aᴵ[i,i] * J)
        zero!(ṽ)
        for m = 1:K
            @← f_ṽ = f(ṽ, t[n] + h * cᴵ[i])
            @. Δṽ = v + h * Aᴵ[i,i] * f_ṽ - ṽ
            ldiv!(Z, Δṽ)
            if norm(Δṽ) < ϵᵣ * norm(ṽ)
                break
            end
            @. ṽ += Δṽ
        end
        @← kᴵ[i] = f(ṽ, t[n] + h * cᴵ[i])
        @← kᴱ[i] = g(ṽ, t[n] + h * cᴱ[i])
    end
    zero!(v)
    for i = 1:s
        @. v += bᴵ[i] * kᴵ[i] + bᴱ[i] * kᴱ[i]
    end
    @. v = u[n] + h * v
    t[n+1] = t[n] + h
end
