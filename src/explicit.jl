"""
    ExplicitRungeKuttaSolver(tableau, p, p̂, h, adaptive) -> RungeKuttaSolver
    ERK(args...; kwargs...) -> RungeKuttaSolver

returns a constructor for an explicit `RungeKuttaSolver` with
- `tableau` : Butcher tableau.
- `p` : order of accuracy.
- `p̂` : order of accuracy of embedded method.
- `h` : step-size.
- `adaptive`: true/false or set of parameters.
"""
mutable struct ExplicitRungeKuttaSolver{tableau_T, p_T, p̂_T, h_T, adaptive_T} <: RungeKuttaSolver
    tableau::tableau_T
    p::p_T
    p̂::p̂_T
    h::h_T
    adaptive::adaptive_T
end

@doc (@doc ExplicitRungeKuttaSolver) ERK(args...; kwargs...) = ExplicitRungeKuttaSolver(args...; kwargs...)

function Base.copy(solver::ExplicitRungeKuttaSolver)
    @↓ tableau, p, p̂, h, adaptive = solver
    return ERK(tableau, p, p̂, h, adaptive)
end

# ---------------------------------- Euler ----------------------------------- #

@doc raw"""
    Euler(; h = 0.0, adaptive = false) -> ExplicitRungeKuttaSolver
    ExplicitEuler(args...; kwargs...) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the <u>explicit</u> <u>1st-order</u> Euler method:
```math
    u_{n+1} = u_n + hf(u_n,~t_n).
```
"""
function Euler(; h = 0.0, adaptive = false)
    tableau = ButcherTableau(
        [0. 0.;
         0. 1.]
    )
    p = 1
    return ERK(tableau, p, ∅, h, adaptive)
end
@doc (@doc Euler) ExplicitEuler(args...; kwargs...) = Euler(args...; kwargs...)

# --------------------------------- Midpoint --------------------------------- #

@doc raw"""
    Midpoint(; h = 0.0, adaptive = false) -> ExplicitRungeKuttaSolver
    ExplicitMidpoint(args...; kwargs...) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the <u>explicit</u> <u>2nd-order</u> Midpoint method:
```math
    u_{n+1} = u_n + k_2
```
with
```math
\begin{align}
    k_1 &= hf(u_n, ~t_n), \\
    k_2 &= hf\textstyle(u_n + \frac{1}{2}k_1, ~t_n + \frac{1}{2}h).
\end{align}
```
"""
function Midpoint(; h = 0.0, adaptive = false)
    tableau = ButcherTableau(
        [ 0.  0.  0.;
         1/2 1/2  0.;
          0.  0.  1.]
    )
    p = 2
    return ERK(tableau, p, ∅, h, adaptive)
end
@doc (@doc Midpoint) ExplicitMidpoint(args...; kwargs...) = Midpoint(args...; kwargs...)

# --------------------------------- Heun 2nd --------------------------------- #

@doc raw"""
    Heun2(; h = 0.0, adaptive = false) -> ExplicitRungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the <u>explicit</u> <u>2nd-order</u> Heun method:
```math
    u_{n+1} = u_n + \textstyle \frac{1}{2} (k_1 + k_2)
```
with
```math
\begin{align}
    k_1 &= hf(u_n, ~t_n), \\
    k_2 &= hf(u_n + k_1, ~t_n + h).
\end{align}
```
"""
function Heun2(; h = 0.0, adaptive = false)
    tableau = ButcherTableau(
        [0.  0.  0.;
         1.  1.  0.;
         0. 1/2 1/2]
    )
    p = 2
    return ERK(tableau, p, ∅, h, adaptive)
end

# ------------------------------- Ralston 2nd -------------------------------- #

@doc raw"""
    Ralston2(; h = 0.0, adaptive = false) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the <u>explicit</u> <u>2nd-order</u> Ralston method:
```math
    u_{n+1} = u_n + \textstyle \frac{1}{4} (k_1 + 3k_2)
```
with
```math
\begin{align}
    k_1 &= hf(u_n, ~t_n), \\
    k_2 &= hf\textstyle(u_n + \frac{2}{3}k_1, ~t_n + \frac{2}{3}h).
\end{align}
```
"""
function Ralston2(; h = 0.0, adaptive = false)
    tableau = ButcherTableau(
        [ 0.  0.  0.;
         2/3 2/3  0.;
          0. 1/4 3/4]
    )
    p = 2
    return ERK(tableau, p, ∅, h, adaptive)
end

# --------------------------------- Heun 3rd --------------------------------- #

@doc raw"""
    Heun3(; h = 0.0, adaptive = false) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the <u>explicit</u> <u>3rd-order</u> Heun method:
```math
    u_{n+1} = u_n + \textstyle \frac{1}{4} (k_1 + 3k_3)
```
with
```math
\begin{align}
    k_1 &= hf(u_n, ~t_n), \\
    k_2 &= hf\textstyle(u_n + \frac{1}{3}k_1, ~t_n + \frac{1}{3}h), \\
    k_3 &= hf\textstyle(u_n + \frac{2}{3}k_1, ~t_n + \frac{2}{3}h).
\end{align}
```
"""
function Heun3(; h = 0.0, adaptive = false)
    tableau = ButcherTableau(
        [ 0.  0.  0.  0.;
         1/3 1/3  0.  0.;
         2/3  0. 2/3  0.;
          0. 1/4  0. 3/4]
    )
    p = 3
    return ERK(tableau, p, ∅, h, adaptive)
end

# -------------------------------- Kutta 3rd --------------------------------- #

@doc raw"""
    Kutta3(; h = 0.0, adaptive = false) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the <u>explicit</u> <u>3rd-order</u> Kutta method:
```math
    u_{n+1} = u_n + \textstyle \frac{1}{6} (k_1 + 3k_2 + k_3)
```
with
```math
\begin{align}
    k_1 &= hf(u_n, ~t_n), \\
    k_2 &= hf\textstyle(u_n + \frac{1}{2}k_1, ~t_n + \frac{1}{2}h), \\
    k_3 &= hf\textstyle(u_n - k_1 + 2k_2, ~t_n + h).
\end{align}
```
"""
function Kutta3(; h = 0.0, adaptive = false)
    tableau = ButcherTableau(
        [ 0.  0.  0.  0.;
         1/2 1/2  0.  0.;
          1. -1.  2.  0.;
          0. 1/6 2/3 1/6]
    )
    p = 3
    return ERK(tableau, p, ∅, h, adaptive)
end

# ------------------------------- Ralston 3rd -------------------------------- #

@doc raw"""
    Ralston3(; h = 0.0, adaptive = false) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the <u>explicit</u> <u>3rd-order</u> Ralston method:
```math
    u_{n+1} = u_n + \textstyle \frac{1}{9} (2k_1 + 3k_2 + 4k_3)
```
with
```math
\begin{align}
    k_1 &= hf(u_n, ~t_n), \\
    k_2 &= hf\textstyle(u_n + \frac{1}{2}k_1, ~t_n + \frac{1}{2}h), \\
    k_3 &= hf\textstyle(u_n + \frac{3}{4}k_1, ~t_n + \frac{3}{4}h).
\end{align}
```
"""
function Ralston3(; h = 0.0, adaptive = false)
    tableau = ButcherTableau(
        [ 0.   0.   0.   0.;
         1/2  1/2   0.   0.;
         3/4   0.  3/4   0.;
          0.  2/9  1/3  4/9]
    )
    p = 3
    return ERK(tableau, p, ∅, h, adaptive)
end

# ----------------------------- Runge-Kutta 4th ------------------------------ #

@doc raw"""
    RK4(; h = 0.0, adaptive = false) -> RungeKuttaSolver

returns an `ExplicitRungeKuttaSolver` for the <u>explicit</u> <u>4th-order</u> Runge-Kutta method:
```math
    u_{n+1} = u_n + \textstyle \frac{1}{6} (k_1 + 2k_2 + 2k_3 + k_4)
```
with
```math
\begin{align}
    k_1 &= hf(u_n, ~t_n), \\
    k_2 &= hf\textstyle(u_n + \frac{1}{2}k_1, ~t_n + \frac{1}{2}h), \\
    k_3 &= hf\textstyle(u_n + \frac{1}{2}k_2, ~t_n + \frac{1}{2}h), \\
    k_4 &= hf(u_n + k_3, ~t_n + h).
\end{align}
```
"""
function RK4(; h = 0.0, adaptive = false)
    tableau = ButcherTableau(
        [ 0.  0.  0.  0.  0.;
         1/2 1/2  0.  0.  0.;
         1/2  0. 1/2  0.  0.;
          1.  0.  0.  1.  0.;
          0. 1/6 1/3 1/3 1/6]
    )
    p = 4
    return ERK(tableau, p, ∅, h, adaptive)
end

# --------------------------- Runge-Kutta-Fehlberg --------------------------- #

@doc raw"""
    RungeKuttaFehlberg(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100)) -> EmbeddedRungeKuttaSolver
    RKF45(args...; kwargs...) -> EmbeddedRungeKuttaSolver

returns a `RungeKuttaSolver` for the Runge-Kutta-Fehlberg method:
```math
\begin{align}
    u_{n+1} &= u_n + \textstyle \frac{25}{216}k_1 + \frac{1408}{2565}k_3 + \frac{2197}{4104}k_4 - \frac{1}{5}k_5, \\
    \tilde{u}_{n+1} &= \tilde{u}_n + \textstyle \frac{16}{135}k_1 + \frac{6656}{12825}k_3 + \frac{28561}{56430}k_4 - \frac{9}{50}k_5 + \frac{2}{55}k_6,
\end{align}
```
where ``\tilde{u}_{n+1}`` is a 5th-order method used to estimate the error in the 4th-order method ``u_{n+1}``, with
```math
\begin{align}
    k_1 &= hf(u_n, ~t_n), \\
    k_2 &= hf\textstyle(u_n + \frac{1}{4}k_1, ~t_n + \frac{1}{4}h), \\
    k_3 &= hf\textstyle(u_n + \frac{3}{32}k_1 + \frac{9}{32}k_2, ~t_n + \frac{3}{8}h), \\
    k_4 &= hf\textstyle(u_n + \frac{1932}{2197}k_1 - \frac{7200}{2197}k_2 + \frac{7296}{2197}k_3, ~t_n + \frac{12}{13}h), \\
    k_5 &= hf\textstyle(u_n + \frac{439}{216}k_1 - 8k_2 + \frac{3680}{513}k_3 - \frac{845}{4104}k_4, ~t_n + h), \\
    k_6 &= hf\textstyle(u_n - \frac{8}{27}k_1 + 2k_2 - \frac{3544}{2565}k_3 + \frac{1859}{4104}k_4 - \frac{11}{40}k_5, ~t_n + \frac{1}{2}h).
\end{align}
```
"""
function RungeKuttaFehlberg(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100))
    tableau = ButcherTableau(
        [   0.          0.          0.          0.          0.          0.          0.;
           1/4         1/4          0.          0.          0.          0.          0.;
           3/8        3/32        9/32          0.          0.          0.          0.;
         12/13   1932/2197  -7200/2197   7296/2197          0.          0.          0.;
            1.     439/216         -8.    3680/513   -845/4104          0.          0.;
           1/2       -8/27          2.  -3544/2565   1859/4104      -11/40          0.;
            0.      25/216          0.   1408/2565   2197/4104        -1/5          0.;
            0.      16/135          0.  6656/12825 28561/56430       -9/50        2/55]
    )
    p = 4
    p̂ = 5
    return ERK(tableau, p, p̂, h, adaptive)
end
@doc (@doc RungeKuttaFehlberg) RKF45(args...; kwargs...) = RungeKuttaFehlberg(args...; kwargs...)

# ------------------------------ Dormand-Price ------------------------------- #

@doc raw"""
    DormandPrince(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100)) -> EmbeddedRungeKuttaSolver
    DOPRI54(args...; kwargs...) -> EmbeddedRungeKuttaSolver

returns a `RungeKuttaSolver` for the Dormand-Prince method:
``\dots``
where ``\tilde{u}_{n+1}`` is a 4th-order method used to estimate the error in the 5th-order method ``u_{n+1}``, with
``\dots``.

"""
function DormandPrince(; h = 0.0, adaptive = AdaptiveParameters(δ = 0.0, ϵ = 1e-5, K = 100))
    tableau = ButcherTableau(
        [  0.            0.            0.            0.            0.            0.            0.            0.;
          1/5           1/5            0.            0.            0.            0.            0.            0.;
         3/10          3/40          9/40            0.            0.            0.            0.            0.;
          4/5         44/45        -56/15          32/9            0.            0.            0.            0.;
          8/9    19372/6561   -25360/2187    64448/6561      -212/729            0.            0.            0.;
           1.     9017/3168       -355/33    46732/5247        49/176   -5103/18656            0.            0.;
           1.        35/384            0.      500/1113       125/192    -2187/6784         11/84            0.;
           0.        35/384            0.      500/1113       125/192    -2187/6784         11/84            0.;
           0.    5179/57600            0.    7571/16695       393/640 -92097/339200      187/2100          1/40]
    )
    p = 5
    p̂ = 4
    return ERK(tableau, p, p̂, h, adaptive)
end
@doc (@doc DormandPrince) DOPRI54(args...; kwargs...) = DormandPrince(args...; kwargs...)
