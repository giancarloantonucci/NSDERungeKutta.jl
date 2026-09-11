# NSDERungeKutta/src/exprk/solvers.jl
#
# Coefficient functions ported VERBATIM from the EXPINT package's scheme files
# (Berland, Skaflestad & Wright, ACM TOMS 33(1), 2007), with Matlab's empty
# entries `[]` rendered as `nothing`. Classical order `p` and stiff order `q`
# are as stated in the EXPINT headers; on semilinear stiff PDEs expect `q`,
# not `p`. `e^z` and `e^{z/2}` are taken from `expphifunctions`' scaling-and-
# squaring chain, NOT reconstructed via `z φ₁(z) + I`: the reconstruction is
# only absolutely accurate and returns 0 or roundoff junk once ‖e^z‖ drops
# below the unit roundoff of ‖zφ₁‖ (strong damping), whereas the chain keeps
# RELATIVE accuracy to the underflow threshold — so stiff modes decay to
# their true e^{hλ} instead of a 1e-16 floor.

"""
    LawsonEuler(; h::Real=0.0) :: ExponentialRungeKuttaSolver

returns an [`ExponentialRungeKuttaSolver`](@ref) for the Lawson–Euler method
(Lawson 1967), classical order 1, stiff order 1.
"""
function LawsonEuler(; h::Real=0.0)
    κ = function (z)
        ez, (φ₁,) = expphifunctions(z, 1)
        Id = one(z)
        U = [Id]
        A = Matrix{Any}(nothing, 1, 1)
        B = Vector{Any}(nothing, 1)
        B[1] = ez
        return U, ez, A, B
    end
    tableau = ExponentialTableau(:LawsonEuler, 1, 1, 1, [0.0], κ)
    return EXPRK(tableau, h)
end

"""
    NorsettEuler(; h::Real=0.0) :: ExponentialRungeKuttaSolver
    ETDEuler(args...; kwargs...) :: ExponentialRungeKuttaSolver
    ExponentialEuler(args...; kwargs...) :: ExponentialRungeKuttaSolver

returns an [`ExponentialRungeKuttaSolver`](@ref) for the Nørsett–Euler method
(also known as ETD-Euler, exponentially fitted Euler, or filtered Euler),
classical order 1, stiff order 1. Exact for constant non-stiff parts at any
step size, since `b₁ = φ₁` integrates a constant integrand exactly.
"""
function NorsettEuler(; h::Real=0.0)
    κ = function (z)
        ez, (φ₁,) = expphifunctions(z, 1)
        Id = one(z)
        U = [Id]
        A = Matrix{Any}(nothing, 1, 1)
        B = Vector{Any}(nothing, 1)
        B[1] = φ₁
        return U, ez, A, B
    end
    tableau = ExponentialTableau(:NorsettEuler, 1, 1, 1, [0.0], κ)
    return EXPRK(tableau, h)
end
@doc (@doc NorsettEuler) ETDEuler(args...; kwargs...) = NorsettEuler(args...; kwargs...)
@doc (@doc NorsettEuler) ExponentialEuler(args...; kwargs...) = NorsettEuler(args...; kwargs...)

"""
    ETD2RK(; h::Real=0.0) :: ExponentialRungeKuttaSolver

returns an [`ExponentialRungeKuttaSolver`](@ref) for the 2-stage ETD2RK
method (Cox & Matthews 2002, eq. 26), classical order 2, stiff order 2.

!!! warning "Deliberate deviation from the EXPINT source"
    EXPINT's shipped `etd2rk.m` (rev. 1.6) sets `b = {φ₁ - 2φ₂, φ₂}`, which
    violates the first order condition `Σᵢ bᵢ(z) = φ₁(z)` — the sum comes to
    `φ₁ - φ₂` — so the shipped scheme is not even consistent (numerically it
    stalls at order 0 with an O(1) bias). The file is a chimera of the
    Hochbruck–Ostermann one-parameter family: `b₁ = φ₁ - 2φ₂` belongs to the
    midpoint variant (`c₂ = ½`, weights `[φ₁ - 2φ₂, 2φ₂]`), while the stages
    it ships (`U₂ = eᶻ`, `a₂₁ = φ₁`, `c₂ = 1`) belong to the endpoint variant
    with weights `[φ₁ - φ₂, φ₂]` — note the file even computes `φ₁(z/2)` and
    never uses it, a leftover of the midpoint form. This port implements the
    endpoint variant with the CORRECT weights, matching Cox & Matthews.
"""
function ETD2RK(; h::Real=0.0)
    κ = function (z)
        ez, (φ₁, φ₂) = expphifunctions(z, 2)
        Id = one(z)
        U = [Id, ez]
        A = Matrix{Any}(nothing, 2, 2)
        A[2,1] = φ₁
        B = Vector{Any}(nothing, 2)
        B[1] = φ₁ - φ₂ # NOT EXPINT's φ₁ - 2φ₂; see the docstring
        B[2] = φ₂
        return U, ez, A, B
    end
    tableau = ExponentialTableau(:ETD2RK, 2, 2, 2, [0.0, 1.0], κ)
    return EXPRK(tableau, h)
end

"""
    ETD3RK(; h::Real=0.0) :: ExponentialRungeKuttaSolver

returns an [`ExponentialRungeKuttaSolver`](@ref) for the 3-stage ETD3RK
method, classical order 3, stiff order 2.
"""
function ETD3RK(; h::Real=0.0)
    κ = function (z)
        ez, (φ₁, φ₂, φ₃) = expphifunctions(z, 3)
        ez2, (φ₁₂,) = expphifunctions(z / 2, 1)
        Id = one(z)
        U = [Id, ez2, ez]
        A = Matrix{Any}(nothing, 3, 3)
        A[2,1] = φ₁₂ / 2
        A[3,1] = -φ₁
        A[3,2] = 2φ₁
        B = Vector{Any}(nothing, 3)
        B[1] = φ₁ - 3φ₂ + 4φ₃
        B[2] = 4φ₂ - 8φ₃
        B[3] = -φ₂ + 4φ₃
        return U, ez, A, B
    end
    tableau = ExponentialTableau(:ETD3RK, 3, 2, 3, [0.0, 0.5, 1.0], κ)
    return EXPRK(tableau, h)
end

"""
    ETD4RK(; h::Real=0.0) :: ExponentialRungeKuttaSolver
    ETDRK4(args...; kwargs...) :: ExponentialRungeKuttaSolver

returns an [`ExponentialRungeKuttaSolver`](@ref) for the classic 4-stage
Cox–Matthews method (Cox & Matthews 2002), classical order 4, stiff order 2
(the well-documented order reduction on stiff semilinear problems; see
Krogstad 2005 and [`HochbruckOstermann4`](@ref) for stiff orders 3 and 4).
"""
function ETD4RK(; h::Real=0.0)
    κ = function (z)
        ez, (φ₁, φ₂, φ₃) = expphifunctions(z, 3)
        ez2, (φ₁₂,) = expphifunctions(z / 2, 1)
        Id = one(z)
        U = [Id, ez2, ez2, ez]
        A = Matrix{Any}(nothing, 4, 4)
        A[2,1] = φ₁₂ / 2
        A[3,2] = φ₁₂ / 2
        A[4,1] = φ₁₂ * φ₁₂ * z / 4
        A[4,3] = φ₁₂
        B = Vector{Any}(nothing, 4)
        B[1] = φ₁ - 3φ₂ + 4φ₃
        B[2] = 2φ₂ - 4φ₃
        B[3] = 2φ₂ - 4φ₃
        B[4] = -φ₂ + 4φ₃
        return U, ez, A, B
    end
    tableau = ExponentialTableau(:ETD4RK, 4, 2, 4, [0.0, 0.5, 0.5, 1.0], κ)
    return EXPRK(tableau, h)
end
@doc (@doc ETD4RK) ETDRK4(args...; kwargs...) = ETD4RK(args...; kwargs...)

"""
    Lawson4(; h::Real=0.0) :: ExponentialRungeKuttaSolver

returns an [`ExponentialRungeKuttaSolver`](@ref) for the classical
Runge–Kutta–Lawson method (Lawson 1967), classical order 4, stiff order 1.
"""
function Lawson4(; h::Real=0.0)
    κ = function (z)
        ez, (φ₁,) = expphifunctions(z, 1)
        ez2, (φ₁₂,) = expphifunctions(z / 2, 1)
        Id = one(z)
        U = [Id, ez2, ez2, ez]
        A = Matrix{Any}(nothing, 4, 4)
        A[2,1] = ez2 / 2
        A[3,2] = Id / 2
        A[4,3] = ez2
        B = Vector{Any}(nothing, 4)
        B[1] = ez / 6
        B[2] = ez2 / 3
        B[3] = ez2 / 3
        B[4] = Id / 6
        return U, ez, A, B
    end
    tableau = ExponentialTableau(:Lawson4, 4, 1, 4, [0.0, 0.5, 0.5, 1.0], κ)
    return EXPRK(tableau, h)
end

"""
    Krogstad(; h::Real=0.0) :: ExponentialRungeKuttaSolver

returns an [`ExponentialRungeKuttaSolver`](@ref) for Krogstad's 4-stage ETD
method (Krogstad 2005), classical order 4, stiff order 3.
"""
function Krogstad(; h::Real=0.0)
    κ = function (z)
        ez, (φ₁, φ₂, φ₃) = expphifunctions(z, 3)
        ez2, (φ₁₂, φ₂₂) = expphifunctions(z / 2, 2)
        Id = one(z)
        U = [Id, ez2, ez2, ez]
        A = Matrix{Any}(nothing, 4, 4)
        A[2,1] = φ₁₂ / 2
        A[3,1] = φ₁₂ / 2 - φ₂₂
        A[3,2] = φ₂₂
        A[4,1] = φ₁ - 2φ₂
        A[4,3] = 2φ₂
        B = Vector{Any}(nothing, 4)
        B[1] = φ₁ - 3φ₂ + 4φ₃
        B[2] = 2φ₂ - 4φ₃
        B[3] = 2φ₂ - 4φ₃
        B[4] = -φ₂ + 4φ₃
        return U, ez, A, B
    end
    tableau = ExponentialTableau(:Krogstad, 4, 3, 4, [0.0, 0.5, 0.5, 1.0], κ)
    return EXPRK(tableau, h)
end

"""
    HochbruckOstermann4(; h::Real=0.0) :: ExponentialRungeKuttaSolver
    HochOst4(args...; kwargs...) :: ExponentialRungeKuttaSolver

returns an [`ExponentialRungeKuttaSolver`](@ref) for the 5-stage exponential
Runge-Kutta method of Hochbruck & Ostermann (2005, p. 19), classical order 4
and STIFF order 4 — the method of choice for stiff semilinear PDEs among the
schemes implemented here.
"""
function HochbruckOstermann4(; h::Real=0.0)
    κ = function (z)
        ez, (φ₁, φ₂, φ₃) = expphifunctions(z, 3)
        ez2, (φ₁₂, φ₂₂, φ₃₂) = expphifunctions(z / 2, 3)
        Id = one(z)
        a₅₂ = φ₂₂ / 2 - φ₃ + φ₂ / 4 - φ₃₂ / 2
        a₅₄ = φ₂₂ / 4 - a₅₂
        U = [Id, ez2, ez2, ez, ez2]
        A = Matrix{Any}(nothing, 5, 5)
        A[2,1] = φ₁₂ / 2
        A[3,1] = φ₁₂ / 2 - φ₂₂
        A[3,2] = φ₂₂
        A[4,1] = φ₁ - 2φ₂
        A[4,2] = φ₂
        A[4,3] = φ₂
        A[5,1] = φ₁₂ / 2 - 2a₅₂ - a₅₄
        A[5,2] = a₅₂
        A[5,3] = a₅₂
        A[5,4] = a₅₄
        B = Vector{Any}(nothing, 5)
        B[1] = φ₁ - 3φ₂ + 4φ₃
        B[4] = -φ₂ + 4φ₃
        B[5] = 4φ₂ - 8φ₃
        return U, ez, A, B
    end
    tableau = ExponentialTableau(:HochbruckOstermann4, 4, 4, 5, [0.0, 0.5, 0.5, 1.0, 0.5], κ)
    return EXPRK(tableau, h)
end
@doc (@doc HochbruckOstermann4) HochOst4(args...; kwargs...) = HochbruckOstermann4(args...; kwargs...)
