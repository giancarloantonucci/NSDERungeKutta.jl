# NSDERungeKutta/src/adaptive.jl

"""
    StepFailure <: SolverFailure

A step cannot advance time or an adaptive step exhausted its nonfinite retry
budget. Invalid user parameters still raise ArgumentError at construction.
"""
struct StepFailure <: SolverFailure
    message::String
end
Base.showerror(io::IO, err::StepFailure) = print(io, "StepFailure: ", err.message)

"""
    AdaptiveParameters <: AbstractAdaptiveParameters

A composite type for the parameters of an adaptive [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
AdaptiveParameters(εₐ::Real=1e-12, εᵣ::Real=1e-5, Mₙ::Integer=100)
```

## Arguments
- `εₐ :: Real` : absolute tolerance
- `εᵣ :: Real` : relative tolerance
- `Mₙ :: Integer` : maximum number of iterations
"""
struct AdaptiveParameters{εₐ_T<:Real,εᵣ_T<:Real,Mₙ_T<:Integer} <: AbstractAdaptiveParameters
    εₐ::εₐ_T
    εᵣ::εᵣ_T
    Mₙ::Mₙ_T
    function AdaptiveParameters(εₐ::Real, εᵣ::Real, Mₙ::Integer)
        isfinite(εₐ) && isfinite(εᵣ) && εₐ >= 0 && εᵣ >= 0 && (εₐ > 0 || εᵣ > 0) ||
            throw(ArgumentError("tolerances must be finite, nonnegative and not both zero"))
        Mₙ >= 1 || throw(ArgumentError("Mₙ must be at least one"))
        return new{typeof(εₐ),typeof(εᵣ),typeof(Mₙ)}(εₐ, εᵣ, Mₙ)
    end
end

function AdaptiveParameters(; εₐ::Real=1e-12, εᵣ::Real=1e-5, Mₙ::Integer=100)
    return AdaptiveParameters(εₐ, εᵣ, Mₙ)
end

#---------------------------------- FUNCTIONS ----------------------------------

"""
    check_adaptive(tableau::AbstractButcherTableau, adaptive) :: Nothing

throws an `ArgumentError` if `adaptive` parameters are given but `tableau`
carries no embedded pair (`d`, `q`). Adaptive stepping needs an embedded
method; failing at construction beats silently stepping at fixed size.
"""
function check_adaptive(
    tableau::AbstractButcherTableau, adaptive::Union{AbstractAdaptiveParameters,Nothing}
)
    if !(adaptive isa Nothing) && (tableau.d isa Nothing || tableau.q isa Nothing)
        throw(
            ArgumentError(
                "`AdaptiveParameters` given, but the tableau has no embedded pair (`d`, `q`). Adaptive stepping is available for the explicit solvers with embedded tableaus (HeunEuler, BogackiShampine, Fehlberg45, DormandPrince54, Verner65, Fehlberg78).",
            ),
        )
    end
    return nothing
end

"""
    adaptivestep!(cache, solution, solver, adaptive::AbstractAdaptiveParameters)

the generic embedded-pair step-size controller. It reads only the stages `k`
and counters from the cache and the embedded weights `d` from the tableau, so
it is family-agnostic: any solver whose tableau carries an embedded pair gets
adaptivity from this one method. (Construction-time [`check_adaptive`](@ref)
guarantees `d` and `q` are present whenever this is reached.) The error is
`RMS(h*Σ(b-d)k / (εₐ + εᵣ*max(abs(u_old),abs(u_new))))`; accept at most one.
A zero component scale accepts only a zero estimated error. The controller uses
safety 0.9 and a step multiplier in [0.2,2]. Finite errors at the retry limit are
accepted with a warning; such steps do not meet the requested tolerance.
"""
function adaptivestep!(
    cache::AbstractRungeKuttaCache,
    solution::AbstractRungeKuttaSolution,
    solver::AbstractRungeKuttaSolver,
    adaptive::AbstractAdaptiveParameters,
    h::Real=solver.stepsize.h,
)
    @↓ n, m, v, k = cache
    @↓ u, t = solution
    @↓ s, b, p, d, q = solver.tableau
    hs = solution.stepsizes
    @↓ εₐ, εᵣ, Mₙ = adaptive

    save_stepsizes = !isnothing(hs)

    # Local step error h Σ(b-d)k, scaled componentwise by both endpoint states.
    # The RMS reduction makes the tolerance invariant under repeated components.
    zero!(v)
    for i in 1:s
        @. v += (b[i] - d[i]) * k[i]
    end
    err2 = zero(real(h * v[firstindex(v)]))
    for j in eachindex(v)
        scale = εₐ + εᵣ * max(abs(u[n][j]), abs(u[n + 1][j]))
        if iszero(scale)
            if !iszero(v[j])
                err2 = oftype(err2, Inf)
                break
            end
        else
            err2 += abs2(h * v[j] / scale)
        end
    end
    err = sqrt(err2 / length(v))
    (isnan(err) || !all(isfinite, u[n + 1])) && (err = oftype(err, Inf))
    r = iszero(err) ? 2.0 : clamp(0.9 * err^(-1 / (min(p, q) + 1)), 0.2, 2.0)
    hused = h # the step the stages were just computed with — what the record must show
    h *= r    # h is from here on the PROPOSAL for the next attempt

    # Check the retry budget before resetting the attempt counter.
    maxed = m ≥ Mₙ && err > 1
    maxed &&
        !isfinite(err) &&
        throw(StepFailure("nonfinite adaptive error estimate exhausted its retry budget"))
    if maxed
        @warn "Maximum number of step-size reductions ($Mₙ) reached at `t = $(t[n])`; accepting the step."
    end

    if err <= 1 || maxed
        if save_stepsizes
            # Record the attempted step, not the proposed next step.
            push!(hs.accepted, hused)
            push!(hs.rejected, eltype(hs.accepted)[]) # open the next step's retry bin
        end
        m = 1
        n += 1
    else
        isfinite(h) && h > 0 && t[n] + typeof(t[n])(h) > t[n] ||
            throw(StepFailure("adaptive step cannot advance time at t = $(t[n])"))
        if save_stepsizes
            push!(hs.rejected[end], hused)
        end
        m += 1
    end

    @↑ cache = n, m
    return h
end
