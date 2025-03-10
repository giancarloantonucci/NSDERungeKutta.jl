function step!(cache::ImplicitExplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::SplitRightHandSide{𝑁, 𝑀}, solver::ImplicitExplicitRungeKuttaSolver) where {𝑁<:NonlinearRightHandSide, 𝑀<:NonlinearRightHandSide}
    @↓ n, e, v, Uᵢ, kᴵ, kᴱ, J = cache
    @↓ u, t = solution
    @↓ fₛ, fₙₛ = rhs
    @↓ Df! = fₛ
    @↓ implicitableau, explicitableau, stepsize, newton = solver
    @↓ Aᴵ ← A, bᴵ ← b, cᴵ ← c, s = implicitableau
    @↓ Aᴱ ← A, bᴱ ← b, cᴱ ← c = explicitableau
    @↓ h = stepsize
    @↓ εᵣ, Mₙ = newton

    # Stages:
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

        # Simplified Newton:
        l = 1
        zero!(Uᵢ)
        ΔUᵢ = kᴱ[i] # to avoid allocs
        # Fᵢ' = I - h * Aᴵ[i,i] * L
        M = factorize(I - h * Aᴵ[i,i] * J)
        for l = 1:Mₙ
            # kᴵ[i] = fₛ(t[n] + h * cᴵ[i], Uᵢ)
            fₛ(kᴵ[i], Uᵢ, t[n] + h * cᴵ[i])
            # Fᵢ = Eᵢ + h * Aᴵ[i,i] * kᴵ[i] - Uᵢ
            @. ΔUᵢ = v + h * Aᴵ[i,i] * kᴵ[i] - Uᵢ
            # ΔUᵢ = Fᵢ' \ Fᵢ
            ldiv!(M, ΔUᵢ)
            # Uᵢ += ΔUᵢ
            @. Uᵢ += ΔUᵢ
            if norm(ΔUᵢ) < εᵣ * norm(Uᵢ)
                break
            end
        end

        # kᴱ[i] = fₙₛ(t[n] + h * cᴱ[i], Uᵢ)
        fₙₛ(kᴱ[i], Uᵢ, t[n] + h * cᴱ[i])
    end

    # Step:
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
    t[n+1] = kahansum(t[n], h, e)

    return u[n+1], t[n+1]
end

function step!(cache::ImplicitExplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, rhs::SplitRightHandSide{𝐿, 𝑁}, solver::ImplicitExplicitRungeKuttaSolver) where {𝐿<:LinearRightHandSide, 𝑁<:NonlinearRightHandSide}
    @↓ n, v, kᴵ, kᴱ, Uᵢ, J, e = cache
    @↓ u, t = solution
    @↓ fₛ, fₙₛ = rhs
    @↓ L, g! = fₛ
    @↓ implicitableau, explicitableau, stepsize = solver
    @↓ Aᴵ ← A, bᴵ ← b, cᴵ ← c, s = implicitableau
    @↓ Aᴱ ← A, bᴱ ← b, cᴱ ← c = explicitableau
    @↓ h = stepsize

    # Stages:
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
        fₛ(kᴵ[i], Uᵢ, t[n] + h * cᴵ[i])
        # kᴱ[i] = fₙₛ(t[n] + h * cᴱ[i], U)
        fₙₛ(kᴱ[i], Uᵢ, t[n] + h * cᴱ[i])
    end

    # Step:
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
    t[n+1] = kahansum(t[n], h, e)

    return u[n+1], t[n+1]
end
