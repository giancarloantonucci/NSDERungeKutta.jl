function step!(cache::ImplicitExplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::SplitRightHandSide{𝑁, 𝑀}, solver::ImplicitExplicitRungeKuttaSolver) where {𝑁<:NonlinearRightHandSide, 𝑀<:NonlinearRightHandSide}
    @↓ n, v, kᴵ, kᴱ, Uᵢ, J, e = cache
    @↓ u, t = solution
    @↓ stiff, nonstiff = rhs
    @↓ Df! = stiff
    @↓ implicitableau, explicitableau, stepsize, newton = solver
    @↓ Aᴵ ← A, bᴵ ← b, cᴵ ← c, s = implicitableau
    @↓ Aᴱ ← A, bᴱ ← b, cᴱ ← c = explicitableau
    @↓ h = stepsize
    @↓ rtol, nits = newton
    # compute stages
    Df!(J, v, u[n], t[n])
    for i = 1:s
        # Eᵢ = u[n] + h * sum(Aᴵ[i,j] * kᴵ[j] + Aᴱ[i,j] * kᴱ[j] for j = 1:i-1)
        zero!(v)
        for j in 1:i-1
            if Aᴵ[i,j] ≠ 0.0
                @. v += Aᴵ[i,j] * kᴵ[j]
            end
            if Aᴱ[i,j] ≠ 0.0
                @. v += Aᴱ[i,j] * kᴱ[j]
            end
        end
        @. v = u[n] + h * v
        # simplified Newton
        ı = 1
        zero!(Uᵢ)
        ΔUᵢ = kᴱ[i] # avoid allocs
        # Fᵢ' = I - h * Aᴵ[i,i] * L
        M = factorize(I - h * Aᴵ[i,i] * J)
        while (ı == 1) || (norm(ΔUᵢ) > rtol * norm(Uᵢ) && ı < nits)
            # kᴵ[i] = fₛ(t[n] + h * cᴵ[i], Uᵢ)
            stiff(kᴵ[i], Uᵢ, t[n] + h * cᴵ[i])
            # Fᵢ = Eᵢ + h * Aᴵ[i,i] * kᴵ[i] - Uᵢ
            @. ΔUᵢ = v + h * Aᴵ[i,i] * kᴵ[i] - Uᵢ
            # ΔUᵢ = Fᵢ' \ Fᵢ
            ldiv!(M, ΔUᵢ)
            # Uᵢ += ΔUᵢ
            @. Uᵢ += ΔUᵢ
            ı += 1
        end
        # kᴱ[i] = fₙₛ(t[n] + h * cᴱ[i], Uᵢ)
        nonstiff(kᴱ[i], Uᵢ, t[n] + h * cᴱ[i])
    end
    # compute step
    # u[n+1] = u[n] + h * sum(bᴵ[i] * kᴵ[i] + bᴱ[i] * kᴱ[i] for i = 1:s)
    zero!(v)
    for i = 1:s
        if bᴵ[i] ≠ 0.0
            @. v += bᴵ[i] * kᴵ[i]
        end
        if bᴱ[i] ≠ 0.0
            @. v += bᴱ[i] * kᴱ[i]
        end
    end
    @. u[n+1] = u[n] + h * v
    # t[n+1] = t[n] + h
    t[n+1] = kahan_sum(t[n], h, e)
    return u[n+1], t[n+1]
end

function step!(cache::ImplicitExplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::SplitRightHandSide{𝐿, 𝑁}, solver::ImplicitExplicitRungeKuttaSolver) where {𝐿<:LinearRightHandSide, 𝑁<:NonlinearRightHandSide}
    @↓ n, v, kᴵ, kᴱ, Uᵢ, J, e = cache
    @↓ u, t = solution
    @↓ stiff, nonstiff = rhs
    @↓ L, g! = stiff
    @↓ implicitableau, explicitableau, stepsize = solver
    @↓ Aᴵ ← A, bᴵ ← b, cᴵ ← c, s = implicitableau
    @↓ Aᴱ ← A, bᴱ ← b, cᴱ ← c = explicitableau
    @↓ h = stepsize
    # compute stages
    for i = 1:s
        # Eᵢ = u[n] + h * sum(Aᴵ[i,j] * kᴵ[j] + Aᴱ[i,j] * kᴱ[j] for j = 1:i-1)
        zero!(v)
        for j in 1:i-1
            if Aᴵ[i, j] ≠ 0.0
                @. v += Aᴵ[i,j] * kᴵ[j]
            end
            if Aᴱ[i, j] ≠ 0.0
                @. v += Aᴱ[i,j] * kᴱ[j]
            end
        end
        @. v = u[n] + h * v
        # Fᵢ' = I - h * Aᴵ[i,i] * L
        M = factorize(I - h * Aᴵ[i,i] * L)
        # Fᵢ = Eᵢ + h * Aᴵ[i,i] * g(t[n] + h * cᴵ[i])
        g! isa Nothing ? zero!(Uᵢ) : g!(Uᵢ, t[n] + h * cᴵ[i])
        @. Uᵢ = v + h * Aᴵ[i,i] * Uᵢ
        # Uᵢ = Fᵢ' \ Fᵢ
        ldiv!(M, Uᵢ)
        # kᴵ[i] = L * Uᵢ + g(t[n] + h * cᴵ[i])
        stiff(kᴵ[i], Uᵢ, t[n] + h * cᴵ[i])
        # kᴱ[i] = fₙₛ(t[n] + h * cᴱ[i], U)
        nonstiff(kᴱ[i], Uᵢ, t[n] + h * cᴱ[i])
    end
    # compute step
    # u[n+1] = u[n] + h * sum(bᴵ[i] * kᴵ[i] + bᴱ[i] * kᴱ[i] for i = 1:s)
    zero!(v)
    for i = 1:s
        if bᴵ[i] ≠ 0.0
            @. v += bᴵ[i] * kᴵ[i]
        end
        if bᴱ[i] ≠ 0.0
            @. v += bᴱ[i] * kᴱ[i]
        end
    end
    @. u[n+1] = u[n] + h * v
    # t[n+1] = t[n] + h
    t[n+1] = kahan_sum(t[n], h, e)
    return u[n+1], t[n+1]
end
