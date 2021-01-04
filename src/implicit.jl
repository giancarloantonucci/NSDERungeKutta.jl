"""
    ImplicitRungeKuttaSolver(tableau, h, ϵ, K[, adaptive]) -> RungeKuttaSolver
    IRK(args...; kwargs...) -> RungeKuttaSolver

returns a constructor for an implicit `RungeKuttaSolver`.

# Arguments
- `tableau  :: ButcherTableau`     : Butcher tableau.
- `h        :: Real`               : step-size.
- `ϵ        :: Real`               : simplified Newton's relative tolerance.
- `K        :: Integer`            : simplified Newton's maximum number of iterations.
- `adaptive :: AdaptiveParameters` : embedded method's parameters.
"""
mutable struct ImplicitRungeKuttaSolver{tableau_T, h_T, ϵ_T, K_T, adaptive_T} <: RungeKuttaSolver
    tableau::tableau_T
    h::h_T
    ϵ::ϵ_T
    K::K_T
    adaptive::adaptive_T
end

function ImplicitRungeKuttaSolver(tableau, h, ϵ, K)
    adaptive = nothing
    return ImplicitRungeKuttaSolver(tableau, h, ϵ, K, adaptive)
end
@doc (@doc ImplicitRungeKuttaSolver) IRK(args...; kwargs...) = ImplicitRungeKuttaSolver(args...; kwargs...)

function Base.copy(solver::ImplicitRungeKuttaSolver)
    @↓ tableau, h, ϵ, K, adaptive = solver
    return IRK(tableau, h, ϵ, K, adaptive)
end

@doc raw"""
    BackwardEuler(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver
    ImplicitEuler(args...; kwargs...) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the backward Euler method:
```math
\begin{array}{c|c}
    1 & 1 \\
    \hline
    1 & 1
\end{array}
```
"""
function BackwardEuler(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau([
        1. 1.;
        1. 1.;
    ])
    return IRK(tableau, h, ϵ, K)
end
@doc (@doc BackwardEuler) ImplicitEuler(args...; kwargs...) = BackwardEuler(args...; kwargs...)

@doc raw"""
    ImplicitMidpoint(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the implicit midpoint method:
```math
\begin{array}{c|c}
    1/2 & 1/2 \\
    \hline
    2   & 1
\end{array}
```
"""
function ImplicitMidpoint(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau([
        1/2 1/2;
        2.  1. ;
    ])
    return IRK(tableau, h, ϵ, K)
end

@doc raw"""
    CrankNicolson(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the Crank-Nicolson method:
```math
\begin{array}{c|cc}
    0   & 0   & 0   \\
    1   & 1/2 & 1/2 \\
    \hline
    2   & 1/2 & 1/2
\end{array}
```
"""
function CrankNicolson(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau([
        0.  0.  0. ;
        1.  1/2 1/2;
        2.  1/2 1/2;
    ])
    return IRK(tableau, h, ϵ, K)
end

@doc raw"""
    SDIRK3(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 3rd-order SDIRK method:
```math
\begin{array}{c|cc}
    \gamma     & \gamma     & 0          \\
    1-\gamma   & 1-2\gamma  & \gamma     \\
    \hline
    3          & 1/2        & 1/2
\end{array}
```
where ``\gamma = (3 + √3)/6``.
"""
function SDIRK3(; h = 0.0, ϵ = 1e-3, K = 10)
    γ = 1/2 + √3/6
    tableau = ButcherTableau([
        γ    γ    0. ;
        1-γ  1-2γ γ  ;
        3.   1/2  1/2;
    ])
    return IRK(tableau, h, ϵ, K)
end

@doc raw"""
    HammerHollingsworth(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver
    HH4(args...; kwargs...) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 4th-order Hammer and Hollingsworth method:
```math
\begin{array}{c|cc}
    \frac{1}{2} - \frac{\sqrt{3}}{6}  & \frac{1}{4}                      & \frac{1}{4} - \frac{\sqrt{3}}{6} \\
    \frac{1}{2} + \frac{\sqrt{3}}{6}  & \frac{1}{4} + \frac{\sqrt{3}}{6} & \frac{1}{4}                      \\
    \hline
    4                                 & 1/2                              & 1/2
\end{array}
```
"""
function HammerHollingsworth(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau([
        1/2-√3/6  1/4       1/4-√3/6  ;
        1/2+√3/6  1/4+√3/6  1/4       ;
        4.        1/2       1/2
    ])
    return IRK(tableau, h, ϵ, K)
end
@doc (@doc HammerHollingsworth) HH4(args...; kwargs...) = HammerHollingsworth(args...; kwargs...)

@doc raw"""
    LobattoIIIA4(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 4th-order Lobatto IIIA method:
```math
\begin{array}{c|ccc}
    0     & 0     & 0     & 0     \\
    1/2   & 5/24  & 1/3   & -1/24 \\
    1     & 1/6   & 2/3   & 1/6   \\
    \hline
    4     & 1/6   & 2/3   & 1/6
\end{array}
```
"""
function LobattoIIIA4(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau([
        0.   0.   0.   0.   ;
        1/2  5/24 1/3  -1/24;
        1.   1/6  2/3  1/6  ;
        4.   1/6  2/3  1/6  ;
    ])
    return IRK(tableau, h, ϵ, K)
end

@doc raw"""
    RadauIIA5(; h = 0.0, ϵ = 1e-3, K = 10) -> ImplicitRungeKuttaSolver
    
returns an `ImplicitRungeKuttaSolver` for the 5th-order Radau IIA method:
```math
\begin{array}{c|ccc}
    \frac{2}{5} - \frac{\sqrt{6}}{10}          & \frac{11}{45} - \frac{7\sqrt{6}}{360}     & \frac{37}{225} - \frac{169\sqrt{6}}{1800} & -\frac{2}{225} + \frac{\sqrt{6}}{75}      \\
    \frac{2}{5} + \frac{\sqrt{6}}{10}          & \frac{37}{225} + \frac{169\sqrt{6}}{1800} & \frac{11}{45} + \frac{7\sqrt{6}}{360}     & -\frac{2}{225} - \frac{\sqrt{6}}{75}      \\
    1                                          & \frac{4}{9} - \frac{\sqrt{6}}{36}         & \frac{4}{9} + \frac{\sqrt{6}}{36}         & \frac{1}{9}                               \\
    \hline
    5                                          & \frac{4}{9} - \frac{\sqrt{6}}{36}         & \frac{4}{9} + \frac{\sqrt{6}}{36}         & \frac{1}{9}
\end{array}

```
"""
function RadauIIA5(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau([
        2/5-√6/10          11/45-7*√6/360     37/225-169*√6/1800 -2/225+√6/75      ;
        2/5+√6/10          37/225+169*√6/1800 11/45+7*√6/360     -2/225-√6/75      ;
        1.                 4/9-√6/36          4/9+√6/36          1/9               ;
        5.                 4/9-√6/36          4/9+√6/36          1/9
    ])
    return IRK(tableau, h, ϵ, K)
end
