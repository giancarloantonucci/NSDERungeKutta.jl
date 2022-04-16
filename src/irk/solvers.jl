"""
    LobattoIIIC2(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitRungeKuttaSolver

returns an [`ImplicitRungeKuttaSolver`](@ref) for the 2nd-order Lobatto IIIC method.
"""
function LobattoIIIC2(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 2
    tableau = ButcherTableau(float([
        0 1/2 -1/2;
        1 1/2  1/2;
        p 1/2  1/2;
    ]))
    newton = SimplifiedNewtonParameters(rtol=rtol, nits=nits)
    return IRK(tableau, h, newton)
end

"""
    RadauIA3(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitRungeKuttaSolver

returns an [`ImplicitRungeKuttaSolver`](@ref) for the 3rd-order Radau IA method.
"""
function RadauIA3(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 3
    tableau = ButcherTableau(float([
          0  1/4 -1/4;
        2/3  1/4 5/12;
          p  1/4  3/4;
    ]))
    newton = SimplifiedNewtonParameters(rtol=rtol, nits=nits)
    return IRK(tableau, h, newton)
end

"""
    RadauIIA3(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitRungeKuttaSolver

returns an [`ImplicitRungeKuttaSolver`](@ref) for the 3rd-order Radau IIA method.
"""
function RadauIIA3(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 3
    tableau = ButcherTableau(float([
        1/3 5/12 -1/12;
         1   3/4   1/4;
         p   3/4   1/4;
    ]))
    newton = SimplifiedNewtonParameters(rtol=rtol, nits=nits)
    return IRK(tableau, h, newton)
end

"""
    GaußLegendre4(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitRungeKuttaSolver
    GaussLegendre4(args...; kwargs...) :: ImplicitRungeKuttaSolver

returns an [`ImplicitRungeKuttaSolver`](@ref) for the 4th-order Gauß-Legendre method.
"""
function GaußLegendre4(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 4
    tableau = ButcherTableau(float([
        1/2-√3/6      1/4 1/4-√3/6;
        1/2+√3/6 1/4+√3/6      1/4;
               p      1/2      1/2;
    ]))
    newton = SimplifiedNewtonParameters(rtol=rtol, nits=nits)
    return IRK(tableau, h, newton)
end
@doc (@doc GaußLegendre4) GaussLegendre4(args...; kwargs...) = GaußLegendre4(args...; kwargs...)

"""
    LobattoIIIA4(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitRungeKuttaSolver

returns an [`ImplicitRungeKuttaSolver`](@ref) for the 4th-order Lobatto IIIA method.
"""
function LobattoIIIA4(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 4
    tableau = ButcherTableau(float([
          0    0   0     0;
        1/2 5/24 1/3 -1/24;
          1  1/6 2/3   1/6;
          p  1/6 2/3   1/6;
    ]))
    newton = SimplifiedNewtonParameters(rtol=rtol, nits=nits)
    return IRK(tableau, h, newton)
end

"""
    LobattoIIIB4(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitRungeKuttaSolver

returns an [`ImplicitRungeKuttaSolver`](@ref) for the 4th-order Lobatto IIIB method.
"""
function LobattoIIIB4(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 4
    tableau = ButcherTableau(float([
          0 1/6 -1/6   0;
        1/2 1/6  1/3   0;
          1 1/6  5/6   0;
          p 1/6  2/3 1/6;
    ]))
    newton = SimplifiedNewtonParameters(rtol=rtol, nits=nits)
    return IRK(tableau, h, newton)
end

"""
    LobattoIIIC4(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitRungeKuttaSolver

returns an [`ImplicitRungeKuttaSolver`](@ref) for the 4th-order Lobatto IIIC method.
"""
function LobattoIIIC4(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 4
    tableau = ButcherTableau(float([
          0 1/6 -1/3   1/6;
        1/2 1/6 5/12 -1/12;
          1 1/6  2/3   1/6;
          p 1/6  2/3   1/6;
    ]))
    newton = SimplifiedNewtonParameters(rtol=rtol, nits=nits)
    return IRK(tableau, h, newton)
end

"""
    RadauI5(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitRungeKuttaSolver

returns an [`ImplicitRungeKuttaSolver`](@ref) for the 5th-order Radau I method.
"""
function RadauI5(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 5
    tableau = ButcherTableau(float([
                0         0               0               0;
        (6-√6)/10 (9+√6)/75     (24+√6)/120 (168-73*√6)/600;
        (6+√6)/10 (9-√6)/75 (168+73*√6)/600     (24-√6)/120;
                p       1/9      (16+√6)/36      (16-√6)/36;
    ]))
    newton = SimplifiedNewtonParameters(rtol=rtol, nits=nits)
    return IRK(tableau, h, newton)
end

"""
    RadauIA5(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitRungeKuttaSolver

returns an [`ImplicitRungeKuttaSolver`](@ref) for the 5th-order Radau IA method.
"""
function RadauIA5(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 5
    tableau = ButcherTableau(float([
                0 1/9     (-1-√6)/18     (-1+√6)/18;
        (6-√6)/10 1/9  (88+7*√6)/360 (88-43*√6)/360;
        (6+√6)/10 1/9 (88+43*√6)/360  (88-7*√6)/360;
                p 1/9     (16+√6)/36     (16-√6)/36;
    ]))
    newton = SimplifiedNewtonParameters(rtol=rtol, nits=nits)
    return IRK(tableau, h, newton)
end

"""
    RadauII5(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitRungeKuttaSolver

returns an [`ImplicitRungeKuttaSolver`](@ref) for the 5th-order Radau II method.
"""
function RadauII5(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 5
    tableau = ButcherTableau(float([
        (4-√6)/10    (24-√6)/120 (24-11*√6)/120   0;
        (4+√6)/10 (24+11*√6)/120    (24+√6)/120   0;
                1      (6-√6)/12      (6+√6)/12   0;
                p     (16-√6)/36     (16+√6)/36 1/9;
    ]))
    newton = SimplifiedNewtonParameters(rtol=rtol, nits=nits)
    return IRK(tableau, h, newton)
end

"""
    RadauIIA5(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitRungeKuttaSolver

returns an [`ImplicitRungeKuttaSolver`](@ref) for the 5th-order Radau IIA method.
"""
function RadauIIA5(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 5
    tableau = ButcherTableau(float([
        (4-√6)/10     (88-7*√6)/360 (296-169*√6)/1800 (-2+3*√6)/75;
        (4+√6)/10 (296+169*√6)/1800     (88+7*√6)/360 (-2-3*√6)/75;
                1        (16-√6)/36        (16+√6)/36          1/9;
                p        (16-√6)/36        (16+√6)/36          1/9;
    ]))
    newton = SimplifiedNewtonParameters(rtol=rtol, nits=nits)
    return IRK(tableau, h, newton)
end

"""
    GaußLegendre6(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitRungeKuttaSolver
    GaussLegendre6(args...; kwargs...) :: ImplicitRungeKuttaSolver

returns an [`ImplicitRungeKuttaSolver`](@ref) for the 6th-order Gauß–Legendre method.
"""
function GaußLegendre6(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 6
    tableau = ButcherTableau(float([
        1/2-√15/10        5/36 2/9-√15/15 5/36-√15/30;
               1/2 5/36+√15/24        2/9 5/36-√15/24;
        1/2+√15/10 5/36+√15/30 2/9+√15/15        5/36;
                 p        5/18        4/9        5/18;
    ]))
    newton = SimplifiedNewtonParameters(rtol=rtol, nits=nits)
    return IRK(tableau, h, newton)
end
@doc (@doc GaußLegendre6) GaussLegendre6(args...; kwargs...) = GaußLegendre6(args...; kwargs...)
