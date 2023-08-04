"""
    BackwardEuler(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: DiagonallyImplicitRungeKuttaSolver
    ImplicitEuler(args...; kwargs...) :: DiagonallyImplicitRungeKuttaSolver

returns an [`DiagonallyImplicitRungeKuttaSolver`](@ref) for the 1st-order backward Euler method.
"""
function BackwardEuler(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 1
    tableau = ButcherTableau(float([
        1 1;
        p 1;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return DIRK(tableau, h, newton)
end
@doc (@doc BackwardEuler) ImplicitEuler(args...; kwargs...) = BackwardEuler(args...; kwargs...)

"""
    ImplicitMidpoint(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: DiagonallyImplicitRungeKuttaSolver
    GaussLegendre2(args...; kwargs...) :: DiagonallyImplicitRungeKuttaSolver

returns an [`DiagonallyImplicitRungeKuttaSolver`](@ref) for the 2nd-order implicit midpoint method.
"""
function ImplicitMidpoint(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 2
    tableau = ButcherTableau(float([
        1/2 1/2;
          p   1;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return DIRK(tableau, h, newton)
end
@doc (@doc ImplicitMidpoint) GaussLegendre2(args...; kwargs...) = ImplicitMidpoint(args...; kwargs...)

"""
    SDIRK2(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: DiagonallyImplicitRungeKuttaSolver

returns an [`DiagonallyImplicitRungeKuttaSolver`](@ref) for the 2nd-order SDIRK method.
"""
function SDIRK2(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 2
    tableau = ButcherTableau(float([
        1   1   0;
        0  -1   1;
        p 1/2 1/2;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return DIRK(tableau, h, newton)
end

"""
    LobattoIII2(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: DiagonallyImplicitRungeKuttaSolver

returns an [`DiagonallyImplicitRungeKuttaSolver`](@ref) for the 2nd-order Lobatto III method.
"""
function LobattoIII2(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 2
    tableau = ButcherTableau(float([
        0   0   0;
        1   1   0;
        p 1/2 1/2;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return DIRK(tableau, h, newton)
end

"""
    CrankNicolson(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: DiagonallyImplicitRungeKuttaSolver
    LobattoIIIA2(args...; kwargs...) :: DiagonallyImplicitRungeKuttaSolver

returns an [`DiagonallyImplicitRungeKuttaSolver`](@ref) for the 2nd-order Crank-Nicolson method.
"""
function CrankNicolson(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 2
    tableau = ButcherTableau(float([
        0   0   0;
        1 1/2 1/2;
        p 1/2 1/2;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return DIRK(tableau, h, newton)
end
@doc (@doc CrankNicolson) LobattoIIIA2(args...; kwargs...) = CrankNicolson(args...; kwargs...)

"""
    SDIRK3(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: DiagonallyImplicitRungeKuttaSolver

returns an [`DiagonallyImplicitRungeKuttaSolver`](@ref) for the 3rd-order SDIRK method.
"""
function SDIRK3(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 3
    γ = 1/2 + √3/6
    tableau = ButcherTableau(float([
          γ    γ   0;
        1-γ 1-2γ   γ;
          p  1/2 1/2;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return DIRK(tableau, h, newton)
end

"""
    RadauI3(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: DiagonallyImplicitRungeKuttaSolver

returns an [`DiagonallyImplicitRungeKuttaSolver`](@ref) for the 3rd-order Radau I method.
"""
function RadauI3(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 3
    tableau = ButcherTableau(float([
         0   0   0;
       2/3 1/3 1/3;
         p 1/4 3/4;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return DIRK(tableau, h, newton)
end

"""
    RadauII3(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: DiagonallyImplicitRungeKuttaSolver

returns an [`DiagonallyImplicitRungeKuttaSolver`](@ref) for the 3rd-order Radau II method.
"""
function RadauII3(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 3
    tableau = ButcherTableau(float([
       1/3 1/3   0;
         1   1   0;
         p 3/4 1/4;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return DIRK(tableau, h, newton)
end

"""
    SDIRK4(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: DiagonallyImplicitRungeKuttaSolver

returns an [`DiagonallyImplicitRungeKuttaSolver`](@ref) for the 4th-order SDIRK method.
"""
function SDIRK4(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 4
    tableau = ButcherTableau(float([
          1/4      1/4         0      0      0   0;
          3/4      1/2       1/4      0      0   0;
        11/20    17/50     -1/25    1/4      0   0;
          1/2 371/1360 -137/2720 15/544    1/4   0;
            1    25/24    -49/48 125/16 -85/12 1/4;
            p    25/24    -49/48 125/16 -85/12 1/4;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return DIRK(tableau, h, newton)
end

"""
    LobattoIII4(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: DiagonallyImplicitRungeKuttaSolver

returns an [`DiagonallyImplicitRungeKuttaSolver`](@ref) for the 4th-order Lobatto III method.
"""
function LobattoIII4(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    p = 4
    tableau = ButcherTableau(float([
          0   0   0   0;
        1/2 1/4 1/4   0;
          1   0   1   0;
          p 1/6 2/3 1/6;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return DIRK(tableau, h, newton)
end
