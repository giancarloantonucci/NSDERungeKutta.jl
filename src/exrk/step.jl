function step!(cache::ExplicitExponentialRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::AbstractRightHandSideFunction, solver::ExplicitExponentialRungeKuttaSolver)
    @↓ n, v, w, k, ehL, ehcL, e = cache
    @↓ u, t = solution
    @↓ tableau, stepsize = solver
    @↓ Aᵩ, bᵩ, c, s = tableau
    @↓ h = stepsize
    # compute stages
    v = u[n+1] # avoid allocs
    for i = 1:s
        zero!(v)
        for j = 1:i-1
            if Aᵩ[i,j] ≠ 0.0
                @. v += Aᵩ[i,j] * k[j]
            end
        end
        mul!(w, ehcL, u[n])
        @. v = w + h * v
        rhs(k[i], v, t[n] + h * c[i])
    end
    # compute step
    zero!(v)
    for i = 1:s
        if bᵩ[i] ≠ 0.0
            @. v += bᵩ[i] * k[i]
        end
    end
    mul!(w, ehL, u[n])
    @. v = w + h * v
    t[n+1] = t[n] +ₖ (h, e)
    return u[n+1], t[n+1]
end
