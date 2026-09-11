# NSDERungeKutta/src/adaptive.jl

"""
    AdaptiveParameters <: AbstractAdaptiveParameters

A composite type for the parameters of an adaptive [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
AdaptiveParameters(εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100)
```

## Arguments
- `εₐ :: Real` : absolute tolerance
- `εᵣ :: Real` : relative tolerance
- `Mₙ :: Integer` : maximum number of iterations
"""
struct AdaptiveParameters{εₐ_T<:Real, εᵣ_T<:Real, Mₙ_T<:Integer} <: AbstractAdaptiveParameters
    εₐ :: εₐ_T
    εᵣ :: εᵣ_T
    Mₙ :: Mₙ_T
end

AdaptiveParameters(; εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100) = AdaptiveParameters(εₐ, εᵣ, Mₙ)

#---------------------------------- FUNCTIONS ----------------------------------

"""
    check_adaptive(tableau::AbstractButcherTableau, adaptive) :: Nothing

throws an `ArgumentError` if `adaptive` parameters are given but `tableau`
carries no embedded pair (`d`, `q`). Adaptive stepping needs an embedded
method; failing at construction beats silently stepping at fixed size.
"""
function check_adaptive(tableau::AbstractButcherTableau, adaptive::Union{AbstractAdaptiveParameters,Nothing})
    if !(adaptive isa Nothing) && (tableau.d isa Nothing || tableau.q isa Nothing)
        throw(ArgumentError("`AdaptiveParameters` given, but the tableau has no embedded pair (`d`, `q`). In v1, adaptive stepping is available for the explicit solvers with embedded tableaus (HeunEuler, BogackiShampine, Fehlberg45, DormandPrince54, Verner65, Fehlberg78)."))
    end
    return nothing
end

"""
    adaptivestep!(cache, solution, solver, adaptive::AbstractAdaptiveParameters)

the generic embedded-pair step-size controller. It reads only the stages `k`
and counters from the cache and the embedded weights `d` from the tableau, so
it is family-agnostic: any solver whose tableau carries an embedded pair gets
adaptivity from this one method. (Construction-time [`check_adaptive`](@ref)
guarantees `d` and `q` are present whenever this is reached.)
"""
function adaptivestep!(cache::AbstractRungeKuttaCache, solution::AbstractRungeKuttaSolution, solver::AbstractRungeKuttaSolver, adaptive::AbstractAdaptiveParameters)
    @↓ n, m, v, k = cache
    @↓ u, t = solution
    @↓ s, b, p, d, q = solver.tableau
    @↓ h, hs = solver.stepsize
    @↓ εₐ, εᵣ, Mₙ = adaptive

    save_stepsizes = !isnothing(hs)

    # Error estimate from the embedded difference, δ = ‖Σᵢ (bᵢ - dᵢ) kᵢ‖:
    zero!(v)
    for i = 1:s
        @. v += (b[i] - d[i]) * k[i]
    end
    ε = εₐ + norm(u[n]) * εᵣ
    δ = hairernorm(v)

    # Step-size update (clamped growth); δ = 0 means the estimate is exact to
    # machine precision, so grow at the clamp rather than divide by zero:
    r = iszero(δ) ? 2.0 : max(0.5, min(2.0, (0.35 * ε / δ)^(1 / (min(p, q) + 1))))
    hused = h # the step the stages were just computed with — what the record must show
    h *= r    # h is from here on the PROPOSAL for the next attempt

    if h ≈ zero(h)
        error("Step-size `h` too small at `t = $(t[n])`.")
    end

    # Check the retry budget BEFORE resetting the counter (the old code reset
    # `m` first, so the warning could never fire):
    maxed = m ≥ Mₙ
    if maxed
        @warn "Maximum number of step-size reductions ($Mₙ) reached at `t = $(t[n])`; accepting the step."
    end

    if δ < ε || maxed
        if save_stepsizes
            # The REALISED step: `hs.accepted == diff(solution.t)`. The old
            # code pushed the post-update h, so "accepted" recorded the
            # proposal for the NEXT step — the name lied (found by the
            # thesis suite's AdaptiveStepsizeDemo, which had to fall back to
            # diff(t)). Same fix below for rejections: record the h that
            # FAILED, not the reduced retry it triggered.
            push!(hs.accepted, hused)
            push!(hs.rejected, []) # open the next step's retry bin
        end
        m = 1
        n += 1
    else
        if save_stepsizes
            push!(hs.rejected[end], hused)
        end
        m += 1
    end

    @↑ cache = n, m
    @↑ solver.stepsize = h
    return solution
end
