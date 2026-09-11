# NSDERungeKutta/src/newton.jl

"""
    NewtonParameters <: AbstractNewtonParameters

A composite type for the parameters of simplified Newton.

# Constructors
```julia
NewtonParameters(; εᵣ=1e-8, εₐ=1e-12, Mₙ=10)
```

## Arguments
- `εᵣ :: Real` : relative tolerance
- `εₐ :: Real` : absolute tolerance
- `Mₙ :: Integer` : maximum number of Newton updates

A stage is accepted when its nonlinear RESIDUAL `r` — the amount by which the
current iterate fails the stage equation — satisfies `‖r‖ ≤ εₐ + εᵣ·scale`,
where `scale` is the larger of the iterate's norm and the norm of the
right-hand side it is compared against. The residual is what is checked, not
the size of the last update: a tiny update is also what a stalled or
badly-scaled iteration produces, so on its own it proves nothing. Non-finite
iterates or residuals are never accepted. If `Mₙ` updates pass without
acceptance, the step throws a [`NewtonFailure`](@ref); an unconverged stage is
never used.

Tolerances must be finite and non-negative; `Mₙ ≥ 1`.

The default `εᵣ = 1e-8` is tighter than the `1e-3` of the former
increment-based test. This is a deliberate accuracy choice for a bound that
now applies to the residual of the iterate actually used; its cost relative
to the old test has not been benchmarked.
"""
struct NewtonParameters{εᵣ_T<:Real, εₐ_T<:Real, Mₙ_T<:Integer} <: AbstractNewtonParameters
    εᵣ::εᵣ_T
    εₐ::εₐ_T
    Mₙ::Mₙ_T
    function NewtonParameters(εᵣ::εᵣ_T, εₐ::εₐ_T, Mₙ::Mₙ_T) where {εᵣ_T<:Real, εₐ_T<:Real, Mₙ_T<:Integer}
        isfinite(εᵣ) && εᵣ ≥ 0 && isfinite(εₐ) && εₐ ≥ 0 && Mₙ ≥ 1 ||
            throw(ArgumentError("`NewtonParameters` needs finite εᵣ ≥ 0, finite εₐ ≥ 0 and Mₙ ≥ 1; got εᵣ = $εᵣ, εₐ = $εₐ, Mₙ = $Mₙ."))
        return new{εᵣ_T, εₐ_T, Mₙ_T}(εᵣ, εₐ, Mₙ)
    end
end

NewtonParameters(εᵣ::Real, Mₙ::Integer) = NewtonParameters(εᵣ, 1e-12, Mₙ) # legacy positional form
NewtonParameters(; εᵣ::Real=1e-8, εₐ::Real=1e-12, Mₙ::Integer=10) = NewtonParameters(εᵣ, εₐ, Mₙ)

"""
    NewtonFailure <: Exception

thrown by an implicit Runge-Kutta step when simplified Newton cannot bring the
stage residual within tolerance: either `Mₙ` updates were spent, or the
iteration left the finite range. Carries the time `t` of the step, the
`stage` (0 for a fully coupled IRK system), the number of `updates` taken, the
last `residual` norm and the `tolerance` it failed to meet. Catch it to reject
or shorten the step; do not use the stage.

Adaptive implicit solvers do not yet recover from this on their own: the
exception is raised before the step-size controller sees the step. That is a
known limit of this release, not a promise the controller makes.
"""
struct NewtonFailure{t_T<:Real, r_T<:Real} <: Exception
    t::t_T
    stage::Int
    updates::Int
    residual::r_T
    tolerance::r_T
end

function Base.showerror(io::IO, e::NewtonFailure)
    where_ = e.stage == 0 ? "the coupled stage system" : "stage $(e.stage)"
    print(io, "NewtonFailure: simplified Newton did not converge on ", where_,
          " at t = ", e.t, " after ", e.updates, " updates (residual ",
          e.residual, ", tolerance ", e.tolerance, "). The stage is unusable; ",
          "shorten the step, raise `Mₙ`, or loosen `εᵣ`/`εₐ`.")
end

"""
    newton_tolerance(newton, xnorm, fnorm) :: Real

the residual bound `εₐ + εᵣ·max(xnorm, fnorm)`, or `NaN` if either scale is
not finite (which no residual can satisfy).
"""
function newton_tolerance(newton::AbstractNewtonParameters, xnorm::Real, fnorm::Real)
    isfinite(xnorm) && isfinite(fnorm) || return oftype(float(xnorm), NaN)
    return newton.εₐ + newton.εᵣ * max(xnorm, fnorm)
end

"""
    newton_accept(rnorm, tol) :: Bool

`true` only for a FINITE residual norm within a finite tolerance. `Inf ≤ Inf`
and `NaN` comparisons are both `false` here on purpose.
"""
newton_accept(rnorm::Real, tol::Real) = isfinite(rnorm) && isfinite(tol) && rnorm ≤ tol

"""
    newton_check(converged, rnorm, tol, newton, t, stage, updates)

throws a [`NewtonFailure`](@ref) unless `converged`.
"""
function newton_check(converged::Bool, rnorm::Real, tol::Real, t::Real, stage::Integer, updates::Integer)
    converged && return nothing
    throw(NewtonFailure(t, Int(stage), Int(updates), promote(float(rnorm), float(tol))...))
end

# `norm` of a vector of stage vectors, and `all(isfinite, ·)` for the same.
# Accumulated with `hypot`, never as a sum of squares: squaring a component
# norm of 1e-200 underflows to 0 (a zero residual would then be "accepted"),
# and 1e200 overflows to Inf, while both are ordinary representable values.
stagesnorm(k::AbstractVector{<:AbstractVector}) = mapreduce(norm, hypot, k; init=zero(float(real(eltype(eltype(k))))))
stagesnorm(k::AbstractVector) = norm(k)
stagesfinite(k::AbstractVector{<:AbstractVector}) = all(x -> all(isfinite, x), k)
stagesfinite(k::AbstractVector) = all(isfinite, k)
