"""
    AdditiveRungeKuttaSolver(ℐ, ℰ) -> RungeKuttaSolver
    ARK(args...; kwargs...) -> RungeKuttaSolver

returns a constructor for an implicit-explicit (ARK) `RungeKuttaSolver` with
- `ℐ` : stiff implicit solver.
- `ℰ` : non-stiff explicit solver.
"""
mutable struct AdditiveRungeKuttaSolver{ℐ_T, ℰ_T} <: RungeKuttaSolver
    ℐ::ℐ_T
    ℰ::ℰ_T
end

@doc (@doc AdditiveRungeKuttaSolver) ARK(args...; kwargs...) = AdditiveRungeKuttaSolver(args...; kwargs...)

function Base.copy(solver::AdditiveRungeKuttaSolver)
    @↓ ℐ, ℰ = solver
    return ARK(ℐ, ℰ)
end

# -------------------------------- IMEXEuler --------------------------------- #

@doc raw"""
    IMEXEuler(; h = 0.0, ϵ = 1e-3, K = 10) -> AdditiveRungeKuttaSolver
    AdditiveEuler(args...; kwargs...) -> AdditiveRungeKuttaSolver

returns an `AdditiveRungeKuttaSolver` for the <u>1st-order</u> ARK-Euler method.
"""
function IMEXEuler(; h = 0.0, ϵ = 1e-3, K = 10)
    ℐ = BackwardEuler(h = h, ϵ = ϵ, K = K)
    ℰ = Euler(h = h)
    return ARK(ℐ, ℰ)
end
@doc (@doc IMEXEuler) AdditiveEuler(args...; kwargs...) = IMEXEuler(args...; kwargs...)
