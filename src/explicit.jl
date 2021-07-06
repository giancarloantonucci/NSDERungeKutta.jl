"""
    ExplicitRungeKuttaSolver(tableau, h[, adaptive]) -> RungeKuttaSolver
    ERK(args...; kwargs...) -> RungeKuttaSolver

returns a constructor for an explicit `RungeKuttaSolver`.

# Arguments
- `tableau  :: ButcherTableau`     : Butcher tableau.
- `stepsize :: StepSize`           : step-size.
- `adaptive :: AdaptiveParameters` : embedded method's parameters.
"""
struct ExplicitRungeKuttaSolver{tableau_T, stepsize_T, adaptive_T} <: RungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    adaptive::adaptive_T
end

ExplicitRungeKuttaSolver(tableau, h::Real, adaptive) = ExplicitRungeKuttaSolver(tableau, StepSize(h), adaptive)
ExplicitRungeKuttaSolver(tableau, stepsize) = ExplicitRungeKuttaSolver(tableau, stepsize, nothing)
@doc (@doc ExplicitRungeKuttaSolver) ERK(args...; kwargs...) = ExplicitRungeKuttaSolver(args...; kwargs...)

"""
    Euler(; h = 0.0) -> ExplicitRungeKuttaSolver
    ExplicitEuler(args...; kwargs...) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 1st-order Euler method.
```
"""
function Euler(; h = 0.0)
    tableau = ButcherTableau(typeof(h)[
        0 0;
        1 1;
    ])
    return ERK(tableau, h)
end
@doc (@doc Euler) ExplicitEuler(args...; kwargs...) = Euler(args...; kwargs...)

"""
    Heun2(; h = 0.0) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 2nd-order Heun method.
"""
function Heun2(; h = 0.0)
    tableau = ButcherTableau(typeof(h)[
        0   0   0 ;
        1   1   0 ;
        2  1/2 1/2;
    ])
    return ERK(tableau, h)
end

"""
    Ralston2(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 2nd-order Ralston method.
"""
function Ralston2(; h = 0.0)
    tableau = ButcherTableau(typeof(h)[
         0   0   0 ;
        2/3 2/3  0 ;
         2  1/4 3/4;
    ])
    return ERK(tableau, h)
end

"""
    Midpoint(; h = 0.0) -> ExplicitRungeKuttaSolver
    ExplicitMidpoint(args...; kwargs...) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 2nd-order midpoint method.
"""
function Midpoint(; h = 0.0)
    tableau = ButcherTableau(typeof(h)[
         0   0   0;
        1/2 1/2  0;
         2   0   1;
    ])
    return ERK(tableau, h)
end
@doc (@doc Midpoint) ExplicitMidpoint(args...; kwargs...) = Midpoint(args...; kwargs...)

"""
    Heun3(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 3rd-order Heun method.
"""
function Heun3(; h = 0.0)
    tableau = ButcherTableau(typeof(h)[
         0   0   0   0 ;
        1/3 1/3  0   0 ;
        2/3  0  2/3  0 ;
         3  1/4  0  3/4;
    ])
    return ERK(tableau, h)
end

"""
    Kutta3(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 3rd-order Kutta method.
"""
function Kutta3(; h = 0.0)
    tableau = ButcherTableau(typeof(h)[
         0   0   0   0 ;
        1/2 1/2  0   0 ;
         1  -1   2   0 ;
         3  1/6 2/3 1/6;
    ])
    return ERK(tableau, h)
end

"""
    Ralston3(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 3rd-order Ralston method.
"""
function Ralston3(; h = 0.0)
    tableau = ButcherTableau(typeof(h)[
         0   0   0   0 ;
        1/2 1/2  0   0 ;
        3/4  0  3/4  0 ;
         3  2/9 1/3 4/9;
    ])
    return ERK(tableau, h)
end

"""
    SSPRK3(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 3rd-order Strong-Stability-Preserving Runge-Kutta method.
"""
function SSPRK3(; h = 0.0)
    tableau = ButcherTableau(typeof(h)[
         0   0   0   0 ;
         1   1   0   0 ;
        1/2 1/4 1/4  0 ;
         3  1/6 1/6 2/3;
    ])
    return ERK(tableau, h)
end

# --------------------------- Runge-Kutta, order 4 --------------------------- #

"""
    RK4(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 4th-order Runge-Kutta method.
"""
function RK4(; h = 0.0)
    tableau = ButcherTableau(typeof(h)[
         0   0   0   0   0 ;
        1/2 1/2  0   0   0 ;
        1/2  0  1/2  0   0 ;
         1   0   0   1   0 ;
         4  1/6 1/3 1/3 1/6;
    ])
    return ERK(tableau, h)
end

"""
    Rule38(; h = 0.0) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 4th-order 3/8-rule method.
"""
function Rule38(; h = 0.0)
    tableau = ButcherTableau(typeof(h)[
         0    0    0    0    0 ;
        1/3  1/3   0    0    0 ;
        2/3 -1/3   1    0    0 ;
         1    1   -1    1    0 ;
         4   1/8  3/8  3/8  1/8;
    ])
    return ERK(tableau, h)
end

"""
    HeunEuler(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100)) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 2nd-order Heun-Euler method with 1st-order error estimate.
"""
function HeunEuler(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100))
    tableau = ButcherTableau(typeof(h)[
        0   0   0 ;
        1   1   0 ;
        2  1/2 1/2;
        1   1   0 ;
    ])
    return ERK(tableau, h, adaptive)
end

"""
    Fehlberg45(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100)) -> ExplicitRungeKuttaSolver
    F45(args...; kwargs...) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 4th-order Fehlberg method with 5th-order error estimate.
"""
function Fehlberg45(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100))
    tableau = ButcherTableau(typeof(h)[
          0       0          0          0           0         0    0  ;
         1/4     1/4         0          0           0         0    0  ;
         3/8     3/32       9/32        0           0         0    0  ;
        12/13 1932/2197 -7200/2197  7296/2197       0         0    0  ;
          1    439/216      -8      3680/513    -845/4104     0    0  ;
         1/2    -8/27        2     -3544/2565   1859/4104  -11/40  0  ;
          4     25/216       0      1408/2565   2197/4104   -1/5   0  ;
          5     16/135       0      6656/12825 28561/56430  -9/50 2/55;
    ])
    return ERK(tableau, h, adaptive)
end
@doc (@doc Fehlberg45) F45(args...; kwargs...) = Fehlberg45(args...; kwargs...)

"""
    DormandPrince54(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100)) -> ExplicitRungeKuttaSolver
    DP54(args...; kwargs...) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 5th-order Dormand-Prince method with 4th-order error estimate.
"""
function DormandPrince54(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100))
    tableau = ButcherTableau(typeof(h)[
         0        0            0          0          0          0          0      0  ;
        1/5      1/5           0          0          0          0          0      0  ;
        3/10     3/40         9/40        0          0          0          0      0  ;
        4/5     44/45       -56/15      32/9         0          0          0      0  ;
        8/9  19372/6561  -25360/2187 64448/6561  -212/729       0          0      0  ;
         1    9017/3168    -355/33   46732/5247    49/176  -5103/18656     0      0  ;
         1      35/384         0       500/1113   125/192  -2187/6784    11/84    0  ;
         5      35/384         0       500/1113   125/192  -2187/6784    11/84    0  ;
         4    5179/57600       0      7571/16695  393/640 -92097/339200 187/2100 1/40;
    ])
    return ERK(tableau, h, adaptive)
end
@doc (@doc DormandPrince54) DP54(args...; kwargs...) = DormandPrince54(args...; kwargs...)

"""
    Verner65(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100)) -> ExplicitRungeKuttaSolver
    V65(args...; kwargs...) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the 6th-order Verner method with 5th-order error estimate.
"""
function Verner65(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100))
    tableau = ButcherTableau(typeof(h)[
         0        0          0         0          0          0       0       0        0   ;
        1/6      1/6         0         0          0          0       0       0        0   ;
        4/15     4/75      16/75       0          0          0       0       0        0   ;
        2/3      5/6       -8/3       5/2         0          0       0       0        0   ;
        5/6   -165/64      55/6    -425/64      85/96        0       0       0        0   ;
         1      12/5        -8     4015/612    -11/36      88/255    0       0        0   ;
        1/15 -8263/15000  124/75   -643/680    -81/250   2484/10625  0       0        0   ;
         1    3501/1720  -300/43 297275/52632 -319/2322 24068/84065  0   3850/26703   0   ;
         6       3/40        0      875/2244    23/72     264/1955   0    125/11592 43/616;
         5      13/160       0     2375/5984     5/16      12/85    3/44     0        0   ;
    ])
    return ERK(tableau, h, adaptive)
end
@doc (@doc Verner65) V65(args...; kwargs...) = Verner65(args...; kwargs...)
