"""
step!(cache::RungeKuttaCache, solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver)

computes a step of the `RungeKuttaSolution` of an `InitialValueProblem` using an `ExplicitRungeKuttaSolver`.
"""
function step!(cache::RungeKuttaCache, solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ExplicitRungeKuttaSolver)
    @↓ n, v, k = cache
    @↓ u, t = solution
    @↓ f! = problem.rhs
    @↓ tableau, h = solver
    @↓ s, A, c, b = tableau
    v = view(u, Block(n+1)) # avoid allocs
    # compute stages
    for i = 1:s
        zero!(v)
        for j = 1:i-1
            @. v += A[i,j] * k[Block(j)]
        end
        @. v = u[Block(n)] + h * v
        k_i = view(k, Block(i))
        f!(k_i, v, t[n] + h * c[i])
        t[n+1] = t[n] + h
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
    @↓ n, v, k, Δk, J = cache
    @↓ u, t = solution
    @↓ f!, Df! = problem.rhs
    @↓ tableau, h, ϵ, K = solver
    @↓ s, A, c, b = tableau
    v = view(u, Block(n+1)) # avoid allocs
    Df!(J, v, u[Block(n)], t[n])
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
            Δk_i = view(Δk, Block(i))
            f!(Δk_i, v, t[n] + h * c[i])
            @. Δk_i -= k[Block(i)]
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
