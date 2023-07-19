function step!(cache::DiagonallyImplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::NonlinearRightHandSide, solver::DiagonallyImplicitRungeKuttaSolver)
    @↓ n, v, k, Uᵢ, Δkᵢ, J, e = cache
    @↓ u, t = solution
    @↓ Df! = rhs
    @↓ tableau, stepsize, newton = solver
    @↓ A, b, c, s = tableau
    @↓ h = stepsize
    @↓ rtol, nits = newton
    # compute stages
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
        # simplified Newton
        zero!(k[i])
        # DFᵢ = I - h * A[i,i] * J
        DFᵢ = factorize(I - h * A[i,i] * J)
        for ȷ = 1:nits
            # Uᵢ = Eᵢ + h * A[i,i] * k[i]
            @. Uᵢ = v + h * A[i,i] * k[i]
            # Fᵢ = f(t[n] + h * c[i], Uᵢ) - k[i]
            rhs(Δkᵢ, Uᵢ, t[n] + h * c[i])
            @. Δkᵢ -= k[i]
            # Δkᵢ = DFᵢ \ Fᵢ
            ldiv!(DFᵢ, Δkᵢ)
            # k[i] += Δkᵢ
            @. k[i] += Δkᵢ
            if norm(Δkᵢ) < rtol * norm(k[i])
                break
            end
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
    @. u[n+1] = u[n] + h * v
    # t[n+1] = t[n] + h
    t[n+1] = kahansum(t[n], h, e)
    return u[n+1], t[n+1]
end

function step!(cache::DiagonallyImplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::LinearRightHandSide, solver::DiagonallyImplicitRungeKuttaSolver)
    @↓ n, v, k, e = cache
    @↓ u, t = solution
    @↓ L = rhs
    @↓ tableau, stepsize = solver
    @↓ A, b, c, s = tableau
    @↓ h = stepsize
    # compute stages
    for i = 1:s
        # Eᵢ = u[n] + h * sum(A[i,j] * k[j] for j = 1:i-1)
        zero!(v)
        for j in 1:i-1
            if A[i, j] ≠ 0.0
                @. v += A[i,j] * k[j]
            end
        end
        @. v = u[n] + h * v
        # DFᵢ = I - h * A[i,i] * L
        DFᵢ = factorize(I - h * A[i,i] * L)
        # Fᵢ = L * Eᵢ + g(t[n] + h * c[i])
        rhs(k[i], v, t[n] + h * c[i])
        # k[i] = DFᵢ \ Fᵢ
        ldiv!(DFᵢ, k[i])
    end
    # compute step
    # u[n+1] = u[n] + h * sum(b[i] * k[i] for i = 1:s)
    zero!(v)
    for i = 1:s
        if b[i] ≠ 0.0
            @. v += b[i] * k[i]
        end
    end
    @. u[n+1] = u[n] + h * v
    # t[n+1] = t[n] + h
    t[n+1] = kahansum(t[n], h, e)
    return u[n+1], t[n+1]
end
