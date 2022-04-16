function step!(cache::DiagonallyImplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::NonlinearRightHandSideFunction, solver::DiagonallyImplicitRungeKuttaSolver)
    @↓ n, k, Uᵢ, Δkᵢ, J, e = cache
    @↓ u, t = solution
    @↓ Df! = rhs
    @↓ tableau, stepsize, newton = solver
    @↓ A, b, c, s = tableau
    @↓ h = stepsize
    @↓ rtol, nits = newton
    # compute stages
    v = u[n+1] # avoid allocs
    Df!(J, v, u[n], t[n])
    for i = 1:s
        # Eᵢ = u[n] + h * sum(A[i,j] * k[j] for j = 1:i-1)
        zero!(v)
        for j = 1:i-1
            if A[i,j] ≠ 0.0
                @. v += A[i,j] * k[j]
            end
        end
        @. v = u[n] + h * v
        #simplified Newton
        ı = 1
        zero!(k[i])
        # Fᵢ' = I - h * A[i,i] * J
        M = factorize(I - h * A[i,i] * J)
        while (ı == 1) || (norm(Δkᵢ) > rtol * norm(k[i]) && ı < nits)
            # Uᵢ = Eᵢ + h * A[i,i] * k[i]
            @. Uᵢ = v + h * A[i,i] * k[i]
            # Fᵢ = f(t[n] + h * c[i], Uᵢ) - k[i]
            rhs(Δkᵢ, Uᵢ, t[n] + h * c[i])
            @. Δkᵢ -= k[i]
            # Δkᵢ = Fᵢ' \ Fᵢ
            ldiv!(M, Δkᵢ)
            # k[i] += Δkᵢ
            @. k[i] += Δkᵢ
            ı += 1
        end
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

function step!(cache::DiagonallyImplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::LinearRightHandSideFunction, solver::DiagonallyImplicitRungeKuttaSolver)
    @↓ n, k, e = cache
    @↓ u, t = solution
    @↓ L = rhs
    @↓ tableau, stepsize = solver
    @↓ A, b, c, s = tableau
    @↓ h = stepsize
    # compute stages
    v = u[n+1] # avoid allocs
    for i = 1:s
        # Eᵢ = u[n] + h * sum(A[i,j] * k[j] for j = 1:i-1)
        zero!(v)
        for j in 1:i-1
            if A[i, j] ≠ 0.0
                @. v += A[i,j] * k[j]
            end
        end
        @. v = u[n] + h * v
        # Fᵢ' = I - h * A[i,i] * L
        M = factorize(I - h * A[i,i] * L)
        # Fᵢ = L * Eᵢ + g(t[n] + h * c[i])
        rhs(k[i], v, t[n] + h * c[i])
        # k[i] = Fᵢ' \ Fᵢ
        ldiv!(M, k[i])
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
