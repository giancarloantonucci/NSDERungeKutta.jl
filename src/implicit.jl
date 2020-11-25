"""
    ImplicitRungeKuttaSolver(tableau, p, p̂, h, ϵ, K, adaptive) -> RungeKuttaSolver
    IRK(args...; kwargs...) -> RungeKuttaSolver

returns a constructor for an implicit `RungeKuttaSolver` with
- `tableau` : Butcher tableau.
- `p` : order of accuracy.
- `p̂` : order of accuracy of embedded method.
- `h` : step-size.
- `ϵ` : relative tolerance in simplified Newton's method.
- `K` : maximum number of iterations in simplified Newton's method.
- `adaptive`: true/false or set of parameters.
"""
mutable struct ImplicitRungeKuttaSolver{tableau_T, p_T, p̂_T, h_T, ϵ_T, K_T, adaptive_T} <: RungeKuttaSolver
    tableau::tableau_T
    p::p_T
    p̂::p̂_T
    h::h_T
    ϵ::ϵ_T
    K::K_T
    adaptive::adaptive_T
end

@doc (@doc ImplicitRungeKuttaSolver) IRK(args...; kwargs...) = ImplicitRungeKuttaSolver(args...; kwargs...)

function Base.copy(solver::ImplicitRungeKuttaSolver)
    @↓ tableau, p, p̂, h, ϵ, K, adaptive = solver
    return IRK(tableau, p, p̂, h, ϵ, K, adaptive)
end

# ------------------------------ Backward Euler ------------------------------ #

@doc raw"""
    BackwardEuler(; h = 0.0, ϵ = 1e-3, K = 10, adaptive = false) -> ImplicitRungeKuttaSolver
    ImplicitEuler(args...; kwargs...) -> ImplicitRungeKuttaSolver

returns a `RungeKuttaSolver` for the <u>implicit</u> <u>1st-order</u> Euler method
```math
    u_{n+1} = u_n + hf(u_{n+1}, ~t_{n+1}).
```
"""
function BackwardEuler(; h = 0.0, ϵ = 1e-3, K = 10, adaptive = false)
    tableau = ButcherTableau(
        [1. 1.;
         0. 1.]
    )
    p = 1
    return IRK(tableau, p, ∅, h, ϵ, K, adaptive)
end
@doc (@doc BackwardEuler) ImplicitEuler(args...; kwargs...) = BackwardEuler(args...; kwargs...)

# ---------------------------- Implicit Midpoint ----------------------------- #

@doc raw"""
    ImplicitMidpoint(; h = 0.0, ϵ = 1e-3, K = 10, adaptive = false) -> ImplicitRungeKuttaSolver

returns a `RungeKuttaSolver` for the <u>implicit</u> <u>2nd-order</u> Midpoint method
```math
    u_{n+1} = u_n + hf\textstyle(\frac{1}{2}u_n + \frac{1}{2}u_{n+1}, ~t_n + \frac{1}{2}h\bigr).
```
"""
function ImplicitMidpoint(; h = 0.0, ϵ = 1e-3, K = 10, adaptive = false)
    tableau = ButcherTableau(
        [1/2 1/2;
          0.  1.]
    )
    p = 2
    return IRK(tableau, p, ∅, h, ϵ, K, adaptive)
end

# ------------------------------- Trapezoidal -------------------------------- #

@doc raw"""
    Trapezoidal(; h = 0.0, ϵ = 1e-3, K = 10, adaptive = false) -> ImplicitRungeKuttaSolver
    CrankNicolson(args...; kwargs...) -> ImplicitRungeKuttaSolver

returns a `RungeKuttaSolver` for the <u>implicit</u> <u>2nd-order</u> Trapezoidal method
```math
    u_{n+1} = u_n + \textstyle \frac{1}{2} hf(u_n, ~t_n) + \frac{1}{2} hf(u_{n+1}, ~t_{n+1}).
```
As such, it is a simple average of the `Euler` and `BackwardEuler` schemes.
"""
function Trapezoidal(; h = 0.0, ϵ = 1e-3, K = 10, adaptive = false)
    tableau = ButcherTableau(
        [0.  0.  0.;
         1. 1/2 1/2;
         0. 1/2 1/2]
    )
    p = 2
    return IRK(tableau, p, ∅, h, ϵ, K, adaptive)
end
@doc (@doc Trapezoidal) CrankNicolson(args...; kwargs...) = Trapezoidal(args...; kwargs...)
