function step!(
    cache::ImplicitRungeKuttaCache,
    solution::AbstractRungeKuttaSolution,
    rhs::NonlinearRightHandSide,
    solver::ImplicitRungeKuttaSolver,
    h::Real=solver.stepsize.h,
)
    @↓ n, e, v, V, k, Δk, J = cache
    @↓ u, t = solution
    @↓ Df! = rhs
    @↓ tableau, newton = solver
    @↓ A, b, c, s = tableau
    @↓ Mₙ = newton

    # Stages:
    zero!(k)
    Df!(J, v, u[n], t[n])
    # DF = I - kron(h * A, J)
    DF = factorize(I - h * kron(A, J))
    # Residual-checked simplified Newton on the coupled stage system (see
    # dirk/step.jl for the contract).
    converged = false
    updates = 0
    rnorm = tol = zero(float(h))
    for l in 0:Mₙ
        fnorm = zero(float(h))
        for i in 1:s
            # Uᵢ = u[n] + h * sum(a[i,j] * k[j] for j = 1:s)
            zero!(v)
            for j in 1:s
                if A[i, j] ≠ 0.0
                    @. v += A[i, j] * k[j]
                end
            end
            @. v = u[n] + h * v
            # r[i] = f(t[n] + h * c[i], Uᵢ) - k[i]
            rhs(Δk[i], v, t[n] + h * c[i])
            fnorm = hypot(fnorm, norm(Δk[i])) # stable accumulation, see `stagesnorm`
            @. Δk[i] -= k[i]
        end
        rnorm = stagesnorm(Δk)
        tol = newton_tolerance(newton, stagesnorm(k), fnorm)
        if newton_accept(rnorm, tol)
            converged = true
            break
        end
        (l == Mₙ || !isfinite(rnorm)) && break
        # Δk = DF \ F, through the flat buffer V
        for i in 1:s
            d = length(Δk[i])
            @. V[((i - 1) * d + 1):(i * d)] = Δk[i]
        end
        ldiv!(DF, V)
        for i in 1:s
            d = length(Δk[i])
            @. Δk[i] = V[((i - 1) * d + 1):(i * d)]
        end
        # k += Δk
        for i in eachindex(k)
            @. k[i] += Δk[i]
        end
        updates += 1
        stagesfinite(k) || break
    end
    newton_check(converged, rnorm, tol, t[n], 0, updates) # stage 0: the coupled system

    # Step:
    # u[n+1] = u[n] + h * sum(b[i] * k[i] for i = 1:s)
    zero!(v)
    for i in 1:s
        if b[i] ≠ 0.0
            @. v += b[i] * k[i]
        end
    end
    @. u[n + 1] = u[n] + h * v
    # t[n+1] = t[n] + h
    t[n + 1] = compensated_sum(t[n], h, e)
    return u[n + 1], t[n + 1]
end

function step!(
    cache::ImplicitRungeKuttaCache,
    solution::AbstractRungeKuttaSolution,
    rhs::LinearRightHandSide,
    solver::ImplicitRungeKuttaSolver,
    h::Real=solver.stepsize.h,
)
    @↓ n, v, k, V, e = cache
    @↓ u, t = solution
    @↓ L, g! = rhs
    @↓ tableau = solver
    @↓ A, b, c, s = tableau

    # Stages:
    # DF = I - kron(h * A, L)
    DF = factorize(I - h * kron(A, L))
    # F = kron(ones(s), L * u[n]) + [g(t[n] + h * c[i]) for i = 1:s]
    mul!(v, L, u[n])
    for i in 1:s
        g! isa Nothing ? zero!(k[i]) : g!(k[i], t[n] + h * c[i])
        @. k[i] += v
    end
    # k = DF \ F, through the flat buffer V
    for i in 1:s
        d = length(k[i])
        @. V[((i - 1) * d + 1):(i * d)] = k[i]
    end
    directldiv!(DF, V) # NOT ldiv!: CHOLMOD factors (sparse SPD L) lack it — see utils.jl
    for i in 1:s
        d = length(k[i])
        @. k[i] = V[((i - 1) * d + 1):(i * d)]
    end

    # Step:
    # u[n+1] = u[n] + h * sum(b[i] * k[i] for i = 1:s)
    zero!(v)
    for i in 1:s
        if b[i] ≠ 0.0
            @. v += b[i] * k[i]
        end
    end
    @. u[n + 1] = u[n] + h * v
    # t[n+1] = t[n] + h
    t[n + 1] = compensated_sum(t[n], h, e)

    return u[n + 1], t[n + 1]
end
