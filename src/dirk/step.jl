# NSDERungeKutta/src/dirk/step.jl

function step!(cache::DiagonallyImplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::NonlinearRightHandSide, solver::DiagonallyImplicitRungeKuttaSolver)
    @↓ n, e, v, Uᵢ, Δkᵢ, k, J = cache
    @↓ u, t = solution
    @↓ Df! = rhs
    @↓ tableau, stepsize, newton = solver
    @↓ A, b, c, s = tableau
    @↓ h = stepsize
    @↓ Mₙ = newton

    # Stages:
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

        # Simplified Newton:
        # DFᵢ = I - h * A[i,i] * J
        zero!(k[i])
        DFᵢ = factorize(I - h * A[i,i] * J)
        # Residual-checked simplified Newton: evaluate the stage equation at
        # the current iterate, accept on a small RESIDUAL, otherwise update.
        # `Mₙ` bounds the number of updates; the residual is always evaluated
        # once more after the last update, so an accepted stage has been
        # checked as it stands, and the increment size is never the verdict.
        converged = false
        updates = 0
        rnorm = tol = zero(float(h))
        for l = 0:Mₙ
            # Uᵢ = Eᵢ + h * A[i,i] * k[i]
            @. Uᵢ = v + h * A[i,i] * k[i]
            # rᵢ = f(t[n] + h * c[i], Uᵢ) - k[i]
            rhs(Δkᵢ, Uᵢ, t[n] + h * c[i])
            fnorm = norm(Δkᵢ)
            @. Δkᵢ -= k[i]
            rnorm = norm(Δkᵢ)
            tol = newton_tolerance(newton, norm(k[i]), fnorm)
            if newton_accept(rnorm, tol)
                converged = true
                break
            end
            (l == Mₙ || !isfinite(rnorm)) && break
            # Δkᵢ = DFᵢ \ rᵢ ; k[i] += Δkᵢ
            ldiv!(DFᵢ, Δkᵢ)
            @. k[i] += Δkᵢ
            updates += 1
            all(isfinite, k[i]) || break
        end
        # Never use an unconverged stage: the step is wrong by an unknown
        # amount, and a fixed-step solver has no controller to catch it.
        newton_check(converged, rnorm, tol, t[n], i, updates)
    end

    # Step:
    # u[n+1] = u[n] + h * sum(b[i] * k[i] for i = 1:s)
    zero!(v)
    for i = 1:s
        if b[i] ≠ 0.0
            @. v += b[i] * k[i]
        end
    end
    @. u[n+1] = u[n] + h * v
    # t[n+1] = t[n] + h
    t[n+1] = compensated_sum(t[n], h, e)

    return u[n+1], t[n+1]
end

function step!(cache::DiagonallyImplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::LinearRightHandSide, solver::DiagonallyImplicitRungeKuttaSolver)
    @↓ n, v, Δkᵢ, k, e = cache
    @↓ u, t = solution
    @↓ L = rhs
    @↓ tableau, stepsize = solver
    @↓ A, b, c, s = tableau
    @↓ h = stepsize

    # Stages:
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
        rhs(k[i], Δkᵢ, v, t[n] + h * c[i]) # 4-arg form: `Δkᵢ` is idle here and serves as scratch
        # k[i] = DFᵢ \ Fᵢ
        directldiv!(DFᵢ, k[i]) # NOT ldiv!: CHOLMOD factors (sparse SPD L) lack it — see utils.jl
    end

    # Step:
    # u[n+1] = u[n] + h * sum(b[i] * k[i] for i = 1:s)
    zero!(v)
    for i = 1:s
        if b[i] ≠ 0.0
            @. v += b[i] * k[i]
        end
    end
    @. u[n+1] = u[n] + h * v
    # t[n+1] = t[n] + h
    t[n+1] = compensated_sum(t[n], h, e)

    return u[n+1], t[n+1]
end
