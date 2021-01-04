"""
    ExplicitRungeKuttaSolver(tableau, h[, adaptive]) -> RungeKuttaSolver
    ERK(args...; kwargs...) -> RungeKuttaSolver

returns a constructor for an explicit `RungeKuttaSolver`.

# Arguments
- `tableau  :: ButcherTableau`     : Butcher tableau.
- `h        :: Real`               : step-size.
- `adaptive :: AdaptiveParameters` : embedded method's parameters.
"""
mutable struct ExplicitRungeKuttaSolver{tableau_T, h_T, adaptive_T} <: RungeKuttaSolver
    tableau::tableau_T
    h::h_T
    adaptive::adaptive_T
end

function ExplicitRungeKuttaSolver(tableau, h)
    adaptive = nothing
    return ExplicitRungeKuttaSolver(tableau, h, adaptive)
end
@doc (@doc ExplicitRungeKuttaSolver) ERK(args...; kwargs...) = ExplicitRungeKuttaSolver(args...; kwargs...)

function Base.copy(solver::ExplicitRungeKuttaSolver)
    @↓ tableau, h, adaptive = solver
    return ERK(tableau, h, adaptive)
end

@doc raw"""
    Euler(; h = 0.0) -> ExplicitRungeKuttaSolver
    ExplicitEuler(args...; kwargs...) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the Euler method:
```math
\begin{array}{c|c}
    0 & 0 \\
    \hline
    1 & 1
\end{array}
```
"""
function Euler(; h = 0.0)
    tableau = ButcherTableau([
        0. 0.;
        1. 1.;
    ])
    return ERK(tableau, h)
end
@doc (@doc Euler) ExplicitEuler(args...; kwargs...) = Euler(args...; kwargs...)

@doc raw"""
    Heun2(; h = 0.0) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 2nd-order Heun method:
```math
\begin{array}{c|cc}
    0   & 0   & 0   \\
    1   & 1   & 0   \\
    \hline
    2   & 1/2 & 1/2
\end{array}
```
"""
function Heun2(; h = 0.0)
    tableau = ButcherTableau([
        0.  0.  0. ;
        1.  1.  0. ;
        2.  1/2 1/2;
    ])
    return ERK(tableau, h)
end

@doc raw"""
    Ralston2(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 2nd-order Ralston method:
```math
\begin{array}{c|cc}
    0   & 0   & 0   \\
    2/3 & 2/3 & 0   \\
    \hline
    2   & 1/4 & 3/4
\end{array}
```
"""
function Ralston2(; h = 0.0)
    tableau = ButcherTableau([
        0.  0.  0. ;
        2/3 2/3 0. ;
        2.  1/4 3/4;
    ])
    return ERK(tableau, h)
end

@doc raw"""
    Midpoint(; h = 0.0) -> ExplicitRungeKuttaSolver
    ExplicitMidpoint(args...; kwargs...) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the midpoint method:
```math
\begin{array}{c|cc}
    0   & 0   & 0   \\
    1/2 & 1/2 & 0   \\
    \hline
    2   & 0   & 1
\end{array}
```
"""
function Midpoint(; h = 0.0)
    tableau = ButcherTableau([
        0.  0.  0.;
        1/2 1/2 0.;
        2.  0.  1.;
    ])
    return ERK(tableau, h)
end
@doc (@doc Midpoint) ExplicitMidpoint(args...; kwargs...) = Midpoint(args...; kwargs...)

@doc raw"""
    Heun3(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 3rd-order Heun method:
```math
\begin{array}{c|ccc}
    0   & 0   & 0   & 0   \\
    1/3 & 1/3 & 0   & 0   \\
    2/3 & 0   & 2/3 & 0   \\
    \hline
    3   & 1/4 & 0   & 3/4
\end{array}
```
"""
function Heun3(; h = 0.0)
    tableau = ButcherTableau([
        0.  0.  0.  0. ;
        1/3 1/3 0.  0. ;
        2/3 0.  2/3 0. ;
        3.  1/4 0.  3/4;
    ])
    return ERK(tableau, h)
end

@doc raw"""
    Kutta3(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 3rd-order Kutta method:
```math
\begin{array}{c|ccc}
    0   & 0   & 0   & 0   \\
    1/2 & 1/2 & 0   & 0   \\
    1   & -1  & 2   & 0   \\
    \hline
    3   & 1/6 & 2/3 & 1/6
\end{array}
```
"""
function Kutta3(; h = 0.0)
    tableau = ButcherTableau([
        0.  0.  0.  0. ;
        1/2 1/2 0.  0. ;
        1.  -1. 2.  0. ;
        3.  1/6 2/3 1/6;
    ])
    return ERK(tableau, h)
end

@doc raw"""
    Ralston3(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 3rd-order Ralston method:
```math
\begin{array}{c|ccc}
    0   & 0   & 0   & 0   \\
    1/2 & 1/2 & 0   & 0   \\
    3/4 & 0   & 3/4 & 0   \\
    \hline
    3   & 2/9 & 1/3 & 4/9
\end{array}
```
"""
function Ralston3(; h = 0.0)
    tableau = ButcherTableau([
        0.  0.  0.  0. ;
        1/2 1/2 0.  0. ;
        3/4 0.  3/4 0. ;
        3.  2/9 1/3 4/9;
    ])
    return ERK(tableau, h)
end

@doc raw"""
    SSPRK3(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 3rd-order Strong-Stability-Preserving Runge-Kutta method:
```math
\begin{array}{c|ccc}
    0   & 0   & 0   & 0   \\
    1   & 1   & 0   & 0   \\
    1/2 & 1/4 & 1/4 & 0   \\
    \hline
    3   & 1/6 & 1/6 & 2/3
\end{array}
```
"""
function SSPRK3(; h = 0.0)
    tableau = ButcherTableau([
        0.  0.  0.  0. ;
        1.  1.  0.  0. ;
        1/2 1/4 1/4 0. ;
        3.  1/6 1/6 2/3;
    ])
    return ERK(tableau, h)
end

# --------------------------- Runge-Kutta, order 4 --------------------------- #

@doc raw"""
    RK4(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the classic Runge-Kutta method:
```math
\begin{array}{c|cccc}
    0   & 0   & 0   & 0   & 0   \\
    1/2 & 1/2 & 0   & 0   & 0   \\
    1/2 & 0   & 1/2 & 0   & 0   \\
    1   & 0   & 0   & 1   & 0   \\
    \hline
    4   & 1/6 & 1/3 & 1/3 & 1/6
\end{array}
```
"""
function RK4(; h = 0.0)
    tableau = ButcherTableau([
        0.  0.  0.  0.  0. ;
        1/2 1/2 0.  0.  0. ;
        1/2 0.  1/2 0.  0. ;
        1.  0.  0.  1.  0. ;
        4.  1/6 1/3 1/3 1/6;
    ])
    return ERK(tableau, h)
end

@doc raw"""
    Rule38(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 3/8-rule method:
```math
\begin{array}{c|cccc}
    0    & 0    & 0    & 0    & 0    \\
    1/3  & 1/3  & 0    & 0    & 0    \\
    2/3  & -1/3 & 1    & 0    & 0    \\
    1    & 1    & -1   & 1    & 0    \\
    \hline
    4    & 1/8  & 3/8  & 3/8  & 1/8
\end{array}
```
"""
function Rule38(; h = 0.0)
    tableau = ButcherTableau([
        0.   0.   0.   0.   0. ;
        1/3  1/3  0.   0.   0. ;
        2/3  -1/3 1.   0.   0. ;
        1.   1.   -1.  1.   0. ;
        4.   1/8  3/8  3/8  1/8;
    ])
    return ERK(tableau, h)
end

@doc raw"""
    HeunEuler(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100)) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the Heun-Euler method:
```math
\begin{array}{c|cc}
    0   & 0   & 0   \\
    1   & 1   & 0   \\
    \hline
    2   & 1/2 & 1/2 \\
    1   & 1   & 0
\end{array}
```
"""
function HeunEuler(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100))
    tableau = ButcherTableau([
        0.  0.  0. ;
        1.  1.  0. ;
        2.  1/2 1/2;
        1.  1.  0. ;
    ])
    return ERK(tableau, h, adaptive)
end

@doc raw"""
    Fehlberg45(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100)) -> ExplicitRungeKuttaSolver
    F45(args...; kwargs...) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 4th-order Fehlberg method with 5th-order error estimate:
```math
\begin{array}{r|ccccc}
    0            &             &             &             &             &             &             \\
    1/4          & 1/4         &             &             &             &             &             \\
    3/8          & 3/32        & 9/32        &             &             &             &             \\
    12/13        & 1932/2197   & -7200/2197  & 7296/2197   &             &             &             \\
    1            & 439/216     & -8          & 3680/513    & -845/4104   &             &             \\
    1/2          & -8/27       & 2           & -3544/2565  & 1859/4104   & -11/40      &             \\
    \hline
    4            & 25/216      & 0           & 1408/2565   & 2197/4104   & -1/5        & 0           \\
    5            & 16/135      & 0           & 6656/12825  & 28561/56430 & -9/50       & 2/55
\end{array}
```
"""
function Fehlberg45(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100))
    tableau = ButcherTableau([
        0.          0.          0.          0.          0.          0.          0.  ;
        1/4         1/4         0.          0.          0.          0.          0.  ;
        3/8         3/32        9/32        0.          0.          0.          0.  ;
        12/13       1932/2197   -7200/2197  7296/2197   0.          0.          0.  ;
        1.          439/216     -8.         3680/513    -845/4104   0.          0.  ;
        1/2         -8/27       2.          -3544/2565  1859/4104   -11/40      0.  ;
        4.          25/216      0.          1408/2565   2197/4104   -1/5        0.  ;
        5.          16/135      0.          6656/12825  28561/56430 -9/50       2/55;
    ])
    return ERK(tableau, h, adaptive)
end
@doc (@doc Fehlberg45) F45(args...; kwargs...) = Fehlberg45(args...; kwargs...)

@doc raw"""
    DormandPrince54(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100)) -> ExplicitRungeKuttaSolver
    DP54(args...; kwargs...) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 5th-order Dormand-Prince method with 4th-order error estimate:
```math
\begin{array}{r|ccccc}
    0              &               &               &               &               &               &               &               \\
    1/5            & 1/5           &               &               &               &               &               &               \\
    3/10           & 3/40          & 9/40          &               &               &               &               &               \\
    4/5            & 44/45         & -56/15        & 32/9          &               &               &               &               \\
    8/9            & 19372/6561    & -25360/2187   & 64448/6561    & -212/729      &               &               &               \\
    1              & 9017/3168     & -355/33       & 46732/5247    & 49/176        & -5103/18656   &               &               \\
    1              & 35/384        & 0             & 500/1113      & 125/192       & -2187/6784    & 11/84         &               \\
    \hline
    5              & 35/384        & 0             & 500/1113      & 125/192       & -2187/6784    & 11/84         & 0             \\
    4              & 5179/57600    & 0             & 7571/16695    & 393/640       & -92097/339200 & 187/2100      & 1/40
\end{array}
```
"""
function DormandPrince54(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100))
    tableau = ButcherTableau([
        0.            0.            0.            0.            0.            0.            0.            0.  ;
        1/5           1/5           0.            0.            0.            0.            0.            0.  ;
        3/10          3/40          9/40          0.            0.            0.            0.            0.  ;
        4/5           44/45         -56/15        32/9          0.            0.            0.            0.  ;
        8/9           19372/6561    -25360/2187   64448/6561    -212/729      0.            0.            0.  ;
        1.            9017/3168     -355/33       46732/5247    49/176        -5103/18656   0.            0.  ;
        1.            35/384        0.            500/1113      125/192       -2187/6784    11/84         0.  ;
        5.            35/384        0.            500/1113      125/192       -2187/6784    11/84         0.  ;
        4.            5179/57600    0.            7571/16695    393/640       -92097/339200 187/2100      1/40;
    ])
    return ERK(tableau, h, adaptive)
end
@doc (@doc DormandPrince54) DP54(args...; kwargs...) = DormandPrince54(args...; kwargs...)

@doc raw"""
    Verner65(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100)) -> ExplicitRungeKuttaSolver
    V65(args...; kwargs...) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 6th-order Verner method with 5th-order error estimate:
```math
\begin{array}{r|ccccc}
    0              &               &               &              &               &               &               &               &               \\
    1/6            & 1/6           &               &              &               &               &               &               &               \\
    4/15           & 4/75          & 16/75         &              &               &               &               &               &               \\
    2/3            & 5/6           & -8/3          & 5/2          &               &               &               &               &               \\
    5/6            & -165/64       & 55/6          & -425/64      & 85/96         &               &               &               &               \\
    1              & 12/5          & -8            & 4015/612     & -11/36        & 88/255        &               &               &               \\
    1/15           & -8263/15000   & 124/75        & -643/680     & -81/250       & 2484/10625    & 0             &               &               \\
    1              & 3501/1720     & -300/43       & 297275/52632 & -319/2322     & 24068/84065   & 0             & 3850/26703    &               \\
    \hline
    6              & 3/40          & 0             & 875/2244     & 23/72         & 264/1955      & 0             & 125/11592     & 43/616        \\
    5              & 13/160        & 0             & 2375/5984    & 5/16          & 12/85         & 3/44          & 0             & 0
\end{array}
```
"""
function Verner65(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100))
    tableau = ButcherTableau([
        0.           0.           0.           0.           0.           0.           0.           0.           0.    ;
        1/6          1/6          0.           0.           0.           0.           0.           0.           0.    ;
        4/15         4/75         16/75        0.           0.           0.           0.           0.           0.    ;
        2/3          5/6          -8/3         5/2          0.           0.           0.           0.           0.    ;
        5/6          -165/64      55/6         -425/64      85/96        0.           0.           0.           0.    ;
        1.           12/5         -8.          4015/612     -11/36       88/255       0.           0.           0.    ;
        1/15         -8263/15000  124/75       -643/680     -81/250      2484/10625   0.           0.           0.    ;
        1.           3501/1720    -300/43      297275/52632 -319/2322    24068/84065  0.           3850/26703   0.    ;
        6.           3/40         0.           875/2244     23/72        264/1955     0.           125/11592    43/616;
        5.           13/160       0.           2375/5984    5/16         12/85        3/44         0.           0.    ;
    ])
    return ERK(tableau, h, adaptive)
end
@doc (@doc Verner65) V65(args...; kwargs...) = Verner65(args...; kwargs...)
