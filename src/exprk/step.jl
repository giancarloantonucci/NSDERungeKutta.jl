# NSDERungeKutta/src/exprk/step.jl
#
# EXPINT's `expglm` main loop specialised to one-step methods (r = 1), with
# the stepsize h folded into the 5-argument mul!s instead of pre-scaling the
# stages: Uᵢ = U[i]u[n] + hΣⱼAᵢⱼk[j], k[i] = fₙₛ(Uᵢ, tᵢ) + g(tᵢ), and
# u[n+1] = Vu[n] + hΣᵢBᵢk[i]. All coefficient applications go through mul!,
# so Diagonal (spectral), dense, and 1×1 (scalar-problem) operators share one
# code path with no branching.

function step!(cache::ExponentialRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::SplitRightHandSide{𝐿, 𝑁}, solver::ExponentialRungeKuttaSolver) where {𝐿<:LinearRightHandSide, 𝑁<:NonlinearRightHandSide}
    @↓ n, v, w, e, k, U, V, A, B = cache
    @↓ u, t = solution
    @↓ tableau, stepsize = solver
    @↓ c, s = tableau
    @↓ h = stepsize
    @↓ fₛ, fₙₛ = rhs
    @↓ g! = fₛ

    # Stages:
    for i = 1:s
        # Uᵢ = U[i] * u[n] + h * sum(A[i,j] * k[j] for j = 1:i-1)
        mul!(v, U[i], u[n])
        for j = 1:i-1
            A[i,j] === nothing || mul!(v, A[i,j], k[j], h, true)
        end
        # k[i] = fₙₛ(Uᵢ, tᵢ) + g(tᵢ): the forcing rides with the non-stiff part,
        # exactly as EXPINT treats N(u, t).
        fₙₛ(k[i], v, t[n] + h * c[i])
        if !(g! isa Nothing)
            g!(w, t[n] + h * c[i])
            k[i] .+= w
        end
    end

    # Step:
    # u[n+1] = V * u[n] + h * sum(B[i] * k[i] for i = 1:s)
    mul!(v, V, u[n])
    for i = 1:s
        B[i] === nothing || mul!(v, B[i], k[i], h, true)
    end
    @. u[n+1] = v
    # t[n+1] = t[n] + h
    t[n+1] = compensated_sum(t[n], h, e)

    return u[n+1], t[n+1]
end

# Pure linear problems u' = Lu + g(t): identical structure with fₙₛ ≡ 0, so
# k[i] = g(tᵢ). For g ≡ nothing this collapses to u[n+1] = e^{hL} u[n] — exact
# propagation of the linear flow for EVERY scheme, since V = e^z throughout.
function step!(cache::ExponentialRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::LinearRightHandSide, solver::ExponentialRungeKuttaSolver)
    @↓ n, v, e, k, U, V, A, B = cache
    @↓ u, t = solution
    @↓ tableau, stepsize = solver
    @↓ c, s = tableau
    @↓ h = stepsize
    @↓ g! = rhs

    if g! isa Nothing
        mul!(v, V, u[n])
    else
        for i = 1:s
            g!(k[i], t[n] + h * c[i])
        end
        mul!(v, V, u[n])
        for i = 1:s
            B[i] === nothing || mul!(v, B[i], k[i], h, true)
        end
    end
    @. u[n+1] = v
    t[n+1] = compensated_sum(t[n], h, e)

    return u[n+1], t[n+1]
end
