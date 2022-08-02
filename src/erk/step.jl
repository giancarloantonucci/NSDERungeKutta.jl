function step!(cache::ExplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::AbstractRightHandSide, solver::ExplicitRungeKuttaSolver)
    @↓ n, k, e = cache
    @↓ u, t = solution
    @↓ tableau, stepsize = solver
    @↓ A, b, c, s = tableau
    @↓ h = stepsize
    # compute stages
    v = u[n+1] #avoid allocs
    for i = 1:s
        # Uᵢ = u[n] + h * sum(A[i,j] * k[j] for j = 1:i-1)
        zero!(v)
        for j = 1:i-1
            if A[i,j] ≠ 0.0
                @. v += A[i,j] * k[j]
            end
        end
        @. v = u[n] + h * v
        # k[i] = f(t[n] + h * c[i], Uᵢ)
        rhs(k[i], v, t[n] + h * c[i])
    end
    # compute step
    # u[n+1] = u[n] + h * sum(b[i] * k[i] for i = 1:s)
    zero!(v)
    for i = 1:s
        if b[i] ≠ 0.0
            @. v += b[i] * k[i]
        end
    end
    @. v = u[n] + h * v
    # t[n+1] = t[n] + h
    t[n+1] = t[n] +ₖ (h, e)
    return u[n+1], t[n+1]
end
