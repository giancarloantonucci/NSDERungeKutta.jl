@inline function step!(
    cache::ExplicitRungeKuttaCache,
    solution::AbstractRungeKuttaSolution,
    rhs::AbstractRightHandSide,
    solver::ExplicitRungeKuttaSolver,
    h::Real=solver.stepsize.h,
)
    @↓ n, v, w, e, k = cache
    @↓ u, t = solution
    @↓ tableau = solver
    @↓ A, b, c, s = tableau

    # Stages:
    for i in 1:s
        # Uᵢ = u[n] + h * sum(A[i,j] * k[j] for j = 1:i-1)
        zero!(v)
        for j in 1:(i - 1)
            if A[i, j] ≠ 0.0
                @. v += A[i, j] * k[j]
            end
        end
        @. v = u[n] + h * v
        # k[i] = f(t[n] + h * c[i], Uᵢ)
        rhs(k[i], w, v, t[n] + h * c[i]) # 4-arg form: `w` is scratch, so forcing/split RHSs stay allocation-free
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
