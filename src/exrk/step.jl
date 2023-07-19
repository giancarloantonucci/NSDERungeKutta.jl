function step!(cache::ExplicitExponentialRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::AbstractRightHandSide, solver::ExplicitExponentialRungeKuttaSolver)
    @↓ n, v, w, k, e = cache
    @↓ u, t = solution
    @↓ tableau, stepsize = solver
    @↓ Aᵩ, bᵩ, c, s = tableau
    @↓ nonstiff, stiff = rhs
    @↓ L = nonstiff
    @↓ h = stepsize
    # compute stages
    for i = 1:s
        # Uᵢ = exp(c[i] * h * L) * u[n] + h * sum(Aᵩ[i,j](h * L) * k[j] for j = 1:i-1)
        zero!(v)
        for j = 1:i-1
            if Aᵩ[i,j] ≠ 0.0
                v .+= Aᵩ[i,j](h * L) * k[j]
            end
        end
        e_cᵢhL = exp(c[i] * h * L)
        mul!(w, e_cᵢhL, u[n])
        @. v = w + h * v
        # k[i] = f(t[n] + h * c[i], Uᵢ)
        stiff(k[i], v, t[n] + h * c[i])
    end
    # compute step
    # u[n+1] = exp(h * L) * u[n] + h * sum(bᵩ[i](h * L) * k[i] for i = 1:s)
    zero!(v)
    for i = 1:s
        if bᵩ[i] ≠ 0.0
            @. v += bᵩ[i](h * L) * k[i]
        end
    end
    e_hL = exp(h * L)
    mul!(w, e_hL, u[n])
    @. u[n+1] = w + h * v
    # t[n+1] = t[n] + h
    t[n+1] = kahansum(t[n], h, e)
    return u[n+1], t[n+1]
end
