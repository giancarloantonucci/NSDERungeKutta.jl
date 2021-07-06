"""
    ImplicitRungeKuttaSolver(tableau, h, ϵ, K[, adaptive]) -> RungeKuttaSolver
    IRK(args...; kwargs...) -> RungeKuttaSolver

returns a constructor for an implicit `RungeKuttaSolver`.

# Arguments
- `tableau  :: ButcherTableau`      : Butcher tableau.
- `stepsize :: StepSize`            : step-size.
- `ϵ        :: Real`                : simplified Newton's relative tolerance.
- `K        :: Integer`             : simplified Newton's maximum number of iterations.
- `adaptive :: AdaptiveParameters`  : embedded method's parameters.
"""
struct ImplicitRungeKuttaSolver{tableau_T, stepsize_T, ϵ_T, K_T, adaptive_T} <: RungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    ϵ::ϵ_T
    K::K_T
    adaptive::adaptive_T
end

ImplicitRungeKuttaSolver(tableau, h::Real, ϵ, K, adaptive) = ImplicitRungeKuttaSolver(tableau, StepSize(h), ϵ, K, adaptive)
ImplicitRungeKuttaSolver(tableau, stepsize, ϵ, K) = ImplicitRungeKuttaSolver(tableau, stepsize, ϵ, K, nothing)
@doc (@doc ImplicitRungeKuttaSolver) IRK(args...; kwargs...) = ImplicitRungeKuttaSolver(args...; kwargs...)

"""
    BackwardEuler(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver
    ImplicitEuler(args...; kwargs...) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 1st-order backward Euler method.
"""
function BackwardEuler(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        1 1;
        1 1;
    ])
    return IRK(tableau, h, ϵ, K)
end
@doc (@doc BackwardEuler) ImplicitEuler(args...; kwargs...) = BackwardEuler(args...; kwargs...)

"""
    ImplicitMidpoint(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 2nd-order implicit midpoint method.
"""
function ImplicitMidpoint(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        1/2 1/2;
         2   1 ;
    ])
    return IRK(tableau, h, ϵ, K)
end

"""
    CrankNicolson(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 2nd-order Crank-Nicolson method.
"""
function CrankNicolson(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        0   0   0 ;
        1  1/2 1/2;
        2  1/2 1/2;
    ])
    return IRK(tableau, h, ϵ, K)
end

"""
    SDIRK3(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 3rd-order SDIRK method.
"""
function SDIRK3(; h = 0.0, ϵ = 1e-3, K = 10)
    γ = 1/2 + √3/6
    tableau = ButcherTableau(typeof(h)[
         γ    γ   0 ;
        1-γ 1-2γ  γ ;
         3   1/2 1/2;
    ])
    return IRK(tableau, h, ϵ, K)
end

"""
    GaussLegendre4(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver
    GL4(args...; kwargs...) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 4th-order Gauss-Legendre method.
"""
function GaussLegendre4(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        1/2-√3/6   1/4    1/4-√3/6;
        1/2+√3/6 1/4+√3/6   1/4   ;
           4       1/2      1/2   ;
    ])
    return IRK(tableau, h, ϵ, K)
end
@doc (@doc GaussLegendre4) GL4(args...; kwargs...) = GaussLegendre4(args...; kwargs...)

"""
    GaussLegendre6(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver
    GL6(args...; kwargs...) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 6th-order Gauss–Legendre method.
"""
function GaussLegendre6(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        1/2-√15/10     5/36    2/9-√15/15 5/36-√15/30;
           1/2     5/36+√15/24    2/9     5/36-√15/24;
        1/2+√15/10 5/36+√15/30 2/9+√15/15     5/36   ;
            6          5/18       4/9         5/18   ;
    ])
    return IRK(tableau, h, ϵ, K)
end
@doc (@doc GaussLegendre6) GL6(args...; kwargs...) = GaussLegendre6(args...; kwargs...)

"""
    LobattoIIIA4(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 4th-order Lobatto IIIA method.
"""
function LobattoIIIA4(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
         0   0    0    0  ;
        1/2 5/24 1/3 -1/24;
         1  1/6  2/3  1/6 ;
         4  1/6  2/3  1/6 ;
    ])
    return IRK(tableau, h, ϵ, K)
end

"""
    LobattoIIIB2(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 2nd-order Lobatto IIIB method.
"""
function LobattoIIIB2(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        1/2 1/2  0 ;
        1/2 1/2  0 ;
         2  1/2 1/2;
    ])
    return IRK(tableau, h, ϵ, K)
end

"""
    LobattoIIIB4(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 4th-order Lobatto IIIB method.
"""
function LobattoIIIB4(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
         0  1/6 -1/6  0 ;
        1/2 1/6  1/3  0 ;
         1  1/6  5/6  0 ;
         4  1/6  2/3 1/6;
    ])
    return IRK(tableau, h, ϵ, K)
end

"""
    LobattoIIIC2(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 2nd-order Lobatto IIIC method.
"""
function LobattoIIIC2(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        0  1/2 -1/2;
        1  1/2  1/2;
        2  1/2  1/2;
    ])
    return IRK(tableau, h, ϵ, K)
end

"""
    LobattoIIIC4(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 4th-order Lobatto IIIC method.
"""
function LobattoIIIC4(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
         0   1/6 -1/3   1/6 ;
        1/2  1/6  5/12 -1/12;
         1   1/6  2/3   1/6 ;
         4   1/6  2/3   1/6 ;
    ])
    return IRK(tableau, h, ϵ, K)
end

"""
    RadauIA3(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 3rd-order Radau IA method.
"""
function RadauIA3(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
         0  1/4 -1/4 ;
        2/3 1/4  5/12;
         3  1/4  3/4 ;
    ])
    return IRK(tableau, h, ϵ, K)
end

"""
    RadauIA5(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 5th-order Radau IA method.
"""
function RadauIA5(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
            0     1/9   -1/18-√6/18     -1/18+√6/18  ;
        3/5-√6/10 1/9  11/45+7*√6/360 11/45-43*√6/360;
        3/5+√6/10 1/9 11/45+43*√6/360  11/45-7*√6/360;
            5     1/9    4/9+√6/36       4/9-√6/36   ;
    ])
    return IRK(tableau, h, ϵ, K)
end

"""
    RadauIIA3(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 3rd-order Radau IIA method.
"""
function RadauIIA3(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        1/3  5/12 -1/12;
         1   3/4   1/4 ;
         3   3/4   1/4 ;
    ])
    return IRK(tableau, h, ϵ, K)
end

"""
    RadauIIA5(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 5th-order Radau IIA method.
"""
function RadauIIA5(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        2/5-√6/10   11/45-7*√6/360   37/225-169*√6/1800 -2/225+√6/75;
        2/5+√6/10 37/225+169*√6/1800   11/45+7*√6/360   -2/225-√6/75;
           1           4/9-√6/36          4/9+√6/36          1/9    ;
           5           4/9-√6/36          4/9+√6/36          1/9    ;
    ])
    return IRK(tableau, h, ϵ, K)
end
