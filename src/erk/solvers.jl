"""
    Euler(; h::Real=0.0) :: ExplicitRungeKuttaSolver
    ExplicitEuler(args...; kwargs...) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 1st-order Euler method.
```
"""
function Euler(; h::Real=0.0)
    p = 1
    tableau = ButcherTableau(float([
        0 0;
        p 1;
    ]))
    return ERK(tableau, h)
end
@doc (@doc Euler) ExplicitEuler(args...; kwargs...) = Euler(args...; kwargs...)

"""
    Heun2(; h::Real=0.0) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 2nd-order Heun method.
"""
function Heun2(; h::Real=0.0)
    p = 2
    tableau = ButcherTableau(float([
        0   0   0;
        1   1   0;
        p 1/2 1/2;
    ]))
    return ERK(tableau, h)
end

"""
    Midpoint(; h::Real=0.0) :: ExplicitRungeKuttaSolver
    ExplicitMidpoint(args...; kwargs...) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 2nd-order mid-point method.
"""
function Midpoint(; h::Real=0.0)
    p = 2
    tableau = ButcherTableau(float([
          0   0 0;
        1/2 1/2 0;
          p   0 1;
    ]))
    return ERK(tableau, h)
end
@doc (@doc Midpoint) ExplicitMidpoint(args...; kwargs...) = Midpoint(args...; kwargs...)

"""
    Ralston2(; h::Real=0.0) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 2nd-order Ralston method.
"""
function Ralston2(; h::Real=0.0)
    p = 2
    tableau = ButcherTableau(float([
          0   0   0;
        2/3 2/3   0;
          p 1/4 3/4;
    ]))
    return ERK(tableau, h)
end

"""
    Heun3(; h::Real=0.0) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 3rd-order Heun method.
"""
function Heun3(; h::Real=0.0)
    p = 3
    tableau = ButcherTableau(float([
          0   0   0   0;
        1/3 1/3   0   0;
        2/3   0 2/3   0;
          p 1/4   0 3/4;
    ]))
    return ERK(tableau, h)
end

"""
    RungeKutta3(; h::Real=0.0) :: ExplicitRungeKuttaSolver
    RK3(args...; kwargs...)

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 3rd-order Kutta method.
"""
function RungeKutta3(; h::Real=0.0)
    p = 3
    tableau = ButcherTableau(float([
          0   0   0   0;
        1/2 1/2   0   0;
          1  -1   2   0;
          p 1/6 2/3 1/6;
    ]))
    return ERK(tableau, h)
end
@doc (@doc RungeKutta4) RK3(args...; kwargs...) = RungeKutta3(args...; kwargs...)

"""
    Ralston3(; h::Real=0.0) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 3rd-order Ralston method.
"""
function Ralston3(; h::Real=0.0)
    p = 3
    tableau = ButcherTableau(float([
          0   0   0   0;
        1/2 1/2   0   0;
        3/4   0 3/4   0;
          p 2/9 1/3 4/9;
    ]))
    return ERK(tableau, h)
end

"""
    SSPRK3(; h::Real=0.0) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 3rd-order Strong-Stability-Preserving Runge-Kutta method.
"""
function SSPRK3(; h::Real=0.0)
    p = 3
    tableau = ButcherTableau(float([
          0   0   0   0;
          1   1   0   0;
        1/2 1/4 1/4   0;
          p 1/6 1/6 2/3;
    ]))
    return ERK(tableau, h)
end

"""
    Ralston4(; h::Real=0.0) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 4th-order Ralston method.
"""
function Ralston4(; h::Real=0.0)
    p = 4
    tableau = ButcherTableau(float([
                 0          0           0          0          0;
               0.4        0.4           0          0          0;
        0.45573725 0.29697761  0.15875964          0          0;
                 1 0.21810040 -3.05096516 3.83286476          0;
                 p 0.17476028 -0.55148066 1.20553560 0.17118478;
    ]))
    return ERK(tableau, h)
end

"""
    RungeKutta4(; h::Real=0.0) :: ExplicitRungeKuttaSolver
    RK4(args...; kwargs...)

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 4th-order Runge-Kutta method.
"""
function RungeKutta4(; h::Real=0.0)
    p = 4
    tableau = ButcherTableau(float([
          0   0   0   0   0;
        1/2 1/2   0   0   0;
        1/2   0 1/2   0   0;
          1   0   0   1   0;
          p 1/6 1/3 1/3 1/6;
    ]))
    return ERK(tableau, h)
end
@doc (@doc RungeKutta4) RK4(args...; kwargs...) = RungeKutta4(args...; kwargs...)

"""
    Rule38(; h::Real=0.0) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 4th-order 3/8-rule method.
"""
function Rule38(; h::Real=0.0)
    p = 4
    tableau = ButcherTableau(float([
          0    0   0   0   0;
        1/3  1/3   0   0   0;
        2/3 -1/3   1   0   0;
          1    1  -1   1   0;
          p  1/8 3/8 3/8 1/8;
    ]))
    return ERK(tableau, h)
end

"""
    Butcher5(; h::Real=0.0) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 5th-order Butcher method.
"""
function Butcher5(; h::Real=0.0)
    p = 5
    tableau = ButcherTableau(float([
            0    0    0     0     0     0    0;
          1/4  1/4    0     0     0     0    0;
          1/4  1/8  1/8     0     0     0    0;
          1/2    0    0   1/2     0     0    0;
          3/4 3/16 -3/8   3/8  9/16     0    0;
            1 -3/7  8/7   6/7 -12/7   8/7    0;
            p 7/90    0 16/45  2/15 16/45 7/90;
    ]))
    return ERK(tableau, h)
end

"""
    KuttaNystrom5(; h::Real=0.0) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 5th-order Kutta-Nyström method.
"""
function KuttaNystrom5(; h::Real=0.0)
    p = 5
    tableau = ButcherTableau(float([
            0      0     0       0    0      0       0;
          1/3    1/3     0       0    0      0       0;
          2/5   4/25  6/25       0    0      0       0;
            1    1/4    -3    15/4    0      0       0;
          2/3   2/27  10/9  -50/81 8/81      0       0;
          4/5   2/25 12/25    2/15 8/75      0       0;
            p 23/192     0 125/192    0 -27/64 125/192;
    ]))
    return ERK(tableau, h)
end

"""
    Butcher6(; h::Real=0.0) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 6th-order Butcher method.
"""
function Butcher6(; h::Real=0.0)
    p = 6
    tableau = ButcherTableau(float([
            0        0      0      0       0      0     0      0;
          1/3      1/3      0      0       0      0     0      0;
          2/3        0    2/3      0       0      0     0      0;
          1/3     1/12    1/3  -1/12       0      0     0      0;
          5/6    25/48 -55/24  35/48    15/8      0     0      0;
          1/6     3/20 -11/24   -1/8     1/2   1/10     0      0;
            1 -261/260  33/13 43/156 -118/39 32/195 80/39      0;
            p   13/200      0  11/40   11/40   4/25  4/25 13/200;
    ]))
    return ERK(tableau, h)
end

"""
    Butcher7(; h::Real=0.0) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 7th-order Butcher method.
"""
function Butcher7(; h::Real=0.0)
    p = 7
    tableau = ButcherTableau(float([
            0         0   0        0            0               0         0           0       0      0;
          1/6       1/6   0        0            0               0         0           0       0      0;
          1/3         0 1/3        0            0               0         0           0       0      0;
          1/2       1/8   0      3/8            0               0         0           0       0      0;
         2/11  148/1331   0 150/1331     -56/1331               0         0           0       0      0;
          2/3  -404/243   0  -170/27    4024/1701      10648/1701         0           0       0      0;
          6/7 2466/2401   0 1242/343 -19176/16807    -51909/16807 1053/2401           0       0      0;
            0     5/164   0        0       96/539     -1815/20384 -405/2464     49/1144       0      0;
            1   -113/32   0  -195/22         32/7      29403/3584  -729/512   1029/1408   21/16      0;
            p         0   0        0       32/105 1771561/6289920  243/2560 16807/74880 77/1440 11/270;
    ]))
    return ERK(tableau, h)
end

#####
##### Embedded
#####

"""
    HeunEuler(; h::Real=0.0, εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100, save_stepsizes::Bool=false) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 2nd-order Heun-Euler method with 1st-order error estimate.
"""
function HeunEuler(; h::Real=0.0, εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100, save_stepsizes::Bool=false)
    p = 2
    q = 1
    tableau = ButcherTableau(float([
        0    0   0;
        1    1   0;
        p  1/2 1/2;
        q    1   0;
    ]))
    adaptive = AdaptiveParameters(εₐ=εₐ, εᵣ=εᵣ, Mₙ=Mₙ)
    stepsize = StepSize(h; save_stepsizes)
    return ERK(tableau, stepsize, adaptive)
end

"""
    BogackiShampine(; h::Real=0.0, εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100, save_stepsizes::Bool=false) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 3rd-order Bogacki-Shampine method with 2nd-order error estimate.
"""
function BogackiShampine(; h::Real=0.0, εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100, save_stepsizes::Bool=false)
    p = 3
    q = 2
    tableau = ButcherTableau(float([
          0    0   0   0   0;
        1/2  1/2   0   0   0;
        3/4    0 3/4   0   0;
          1  2/9 1/3 4/9   0;
          p  2/9 1/3 4/9   0;
          q 7/24 1/4 1/3 1/8;
    ]))
    adaptive = AdaptiveParameters(εₐ=εₐ, εᵣ=εᵣ, Mₙ=Mₙ)
    stepsize = StepSize(h; save_stepsizes)
    return ERK(tableau, stepsize, adaptive)
end

"""
    Fehlberg45(; h::Real=0.0, εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100, save_stepsizes::Bool=false) :: ExplicitRungeKuttaSolver
    F45(args...; kwargs...) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 4th-order Fehlberg method with 5th-order error estimate.
"""
function Fehlberg45(; h::Real=0.0, εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100, save_stepsizes::Bool=false)
    p = 4
    q = 5
    tableau = ButcherTableau(float([
            0         0          0          0           0      0    0;
          1/4       1/4          0          0           0      0    0;
          3/8      3/32       9/32          0           0      0    0;
        12/13 1932/2197 -7200/2197  7296/2197           0      0    0;
            1   439/216         -8   3680/513   -845/4104      0    0;
          1/2     -8/27          2 -3544/2565   1859/4104 -11/40    0;
            p    25/216          0  1408/2565   2197/4104   -1/5    0;
            q    16/135          0 6656/12825 28561/56430  -9/50 2/55;
    ]))
    adaptive = AdaptiveParameters(εₐ=εₐ, εᵣ=εᵣ, Mₙ=Mₙ)
    stepsize = StepSize(h; save_stepsizes)
    return ERK(tableau, stepsize, adaptive)
end
@doc (@doc Fehlberg45) F45(args...; kwargs...) = Fehlberg45(args...; kwargs...)

"""
    DormandPrince54(; h::Real=0.0, εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100, save_stepsizes::Bool=false) :: ExplicitRungeKuttaSolver
    DP54(args...; kwargs...) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 5th-order Dormand-Prince method with 4th-order error estimate.
"""
function DormandPrince54(; h::Real=0.0, εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100, save_stepsizes::Bool=false)
    p = 5
    q = 4
    tableau = ButcherTableau(float([
           0          0           0         0         0             0        0    0;
         1/5        1/5           0         0         0             0        0    0;
        3/10       3/40        9/40         0         0             0        0    0;
         4/5      44/45      -56/15      32/9         0             0        0    0;
         8/9 19372/6561 -25360/2187 64448/6561 -212/729             0        0    0;
           1  9017/3168     -355/33 46732/5247   49/176   -5103/18656        0    0;
           1     35/384           0   500/1113  125/192    -2187/6784    11/84    0;
           p     35/384           0   500/1113  125/192    -2187/6784    11/84    0;
           q 5179/57600           0 7571/16695  393/640 -92097/339200 187/2100 1/40;
    ]))
    adaptive = AdaptiveParameters(εₐ=εₐ, εᵣ=εᵣ, Mₙ=Mₙ)
    stepsize = StepSize(h; save_stepsizes)
    return ERK(tableau, stepsize, adaptive)
end
@doc (@doc DormandPrince54) DP54(args...; kwargs...) = DormandPrince54(args...; kwargs...)

"""
    Verner65(; h::Real=0.0, εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100, save_stepsizes::Bool=false) :: ExplicitRungeKuttaSolver
    V65(args...; kwargs...) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 6th-order Verner method with 5th-order error estimate.
"""
function Verner65(; h::Real=0.0, εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100, save_stepsizes::Bool=false)
    p = 6
    q = 5
    tableau = ButcherTableau(float([
           0           0       0            0         0           0    0          0      0;
         1/6         1/6       0            0         0           0    0          0      0;
        4/15        4/75   16/75            0         0           0    0          0      0;
         2/3         5/6    -8/3          5/2         0           0    0          0      0;
         5/6     -165/64    55/6      -425/64     85/96           0    0          0      0;
           1        12/5      -8     4015/612    -11/36      88/255    0          0      0;
        1/15 -8263/15000  124/75     -643/680   -81/250  2484/10625    0          0      0;
           1   3501/1720 -300/43 297275/52632 -319/2322 24068/84065    0 3850/26703      0;
           p        3/40       0     875/2244     23/72    264/1955    0  125/11592 43/616;
           q      13/160       0    2375/5984      5/16       12/85 3/44          0      0;
    ]))
    adaptive = AdaptiveParameters(εₐ=εₐ, εᵣ=εᵣ, Mₙ=Mₙ)
    stepsize = StepSize(h; save_stepsizes)
    return ERK(tableau, stepsize, adaptive)
end
@doc (@doc Verner65) V65(args...; kwargs...) = Verner65(args...; kwargs...)

"""
    Fehlberg78(; h::Real=0.0, εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100, save_stepsizes::Bool=false) :: ExplicitRungeKuttaSolver
    F78(args...; kwargs...) :: ExplicitRungeKuttaSolver

returns an [`ExplicitRungeKuttaSolver`](@ref) for the 7th-order Fehlberg method with 8th-order error estimate.
"""
function Fehlberg78(; h::Real=0.0, εₐ::Real=0.0, εᵣ::Real=1e-5, Mₙ::Integer=100, save_stepsizes::Bool=false)
    p = 7
    q = 8
    tableau = ButcherTableau(float([
           0          0    0      0        0         0       0         0     0      0     0      0      0      0;
        2/27       2/27    0      0        0         0       0         0     0      0     0      0      0      0;
         1/9       1/36 1/12      0        0         0       0         0     0      0     0      0      0      0;
         1/6       1/24    0    1/8        0         0       0         0     0      0     0      0      0      0;
        5/12       5/12    0 -25/16    25/16         0       0         0     0      0     0      0      0      0;
         1/2       1/20    0      0      1/4       1/5       0         0     0      0     0      0      0      0;
         5/6    -25/108    0      0  125/108    -65/27  125/54         0     0      0     0      0      0      0;
         1/6     31/300    0      0        0    61/225    -2/9    13/900     0      0     0      0      0      0;
         2/3          2    0      0    -53/6    704/45  -107/9     67/90     3      0     0      0      0      0;
         1/3    -91/108    0      0   23/108  -976/135  311/54    -19/60  17/6  -1/12     0      0      0      0;
           1  2383/4100    0      0 -341/164 4496/1025 -301/82 2133/4100 45/82 45/164 18/41      0      0      0;
           0      3/205    0      0        0         0   -6/41    -3/205 -3/41   3/41  6/41      0      0      0;
           1 -1777/4100    0      0 -341/164 4496/1025 -289/82 2193/4100 51/82 33/164 12/41      0      1      0;
           p     41/840    0      0        0         0  34/105      9/35  9/35  9/280 9/280 41/840      0      0;
           q          0    0      0        0         0  34/105      9/35  9/35  9/280 9/280      0 41/840 41/840;
    ]))
    adaptive = AdaptiveParameters(εₐ=εₐ, εᵣ=εᵣ, Mₙ=Mₙ)
    stepsize = StepSize(h; save_stepsizes)
    return ERK(tableau, stepsize, adaptive)
end
@doc (@doc Fehlberg78) F78(args...; kwargs...) = Fehlberg78(args...; kwargs...)
