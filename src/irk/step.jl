function step!(cache::ImplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::NonlinearRightHandSide, solver::ImplicitRungeKuttaSolver)
    @↓ n, k, Δk, V, J, e = cache
    @↓ u, t = solution
    @↓ Df! = rhs
    @↓ tableau, newton, stepsize = solver
    @↓ A, b, c, s = tableau
    @↓ h = stepsize
    @↓ rtol, nits = newton
    # compute stages
    ı = 1
    zero!(k)
    v = u[n+1] # avoid allocs
    Df!(J, v, u[n], t[n])
    # F' = I - kron(h * A, J)
    M = factorize(I - h * kron(A, J))
    while (ı == 1) || (norm(Δk) > rtol * norm(k) && ı < nits)
        for i = 1:s
            # Uᵢ = u[n] + h * sum(a[i,j] * k[j] for j = 1:s)
            zero!(v)
            for j in 1:s
                if A[i,j] ≠ 0.0
                    @. v += A[i,j] * k[j]
                end
            end
            @. v = u[n] + h * v
            # F[i] = f(t[n] + h * c[i], Uᵢ) - k[i]
            rhs(Δk[i], v, t[n] + h * c[i])
            @. Δk[i] -= k[i]
        end
        # Δk = F' \ F
        # TO-DO: ldiv!(M, Δk)
        for i = 1:s
            d = length(Δk[i])
            @. V[(i-1)*d+1:i*d] = Δk[i]
        end
        ldiv!(M, V)
        for i = 1:s
            d = length(Δk[i])
            @. Δk[i] = V[(i-1)*d+1:i*d]
        end
        # k += Δk
        # TO-DO: @. k += Δk
        for i in eachindex(k)
            @. k[i] += Δk[i]
        end
        ı += 1
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

function step!(cache::ImplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::LinearRightHandSide, solver::ImplicitRungeKuttaSolver)
    @↓ n, k, V, J, e = cache
    @↓ u, t = solution
    @↓ L, g! = rhs
    @↓ tableau, stepsize = solver
    @↓ A, b, c, s = tableau
    @↓ h = stepsize
    # compute stages
    v = u[n+1] # avoid allocs
    # F' = I - kron(h * A, L)
    M = factorize(I - h * kron(A, L))
    # F = kron(ones(s), L * u[n]) + [g(t[n] + h * c[i]) for i = 1:s]
    mul!(v, L, u[n])
    for i = 1:s
        g! isa Nothing ? zero!(k[i]) : g!(k[i], t[n] + h * c[i])
        @. k[i] += v
    end
    # k = F' \ F
    # TO-DO: ldiv!(M, k)
    for i = 1:s
        d = length(k[i])
        @. V[(i-1)*d+1:i*d] = k[i]
    end
    ldiv!(M, V)
    for i = 1:s
        d = length(k[i])
        @. k[i] = V[(i-1)*d+1:i*d]
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
