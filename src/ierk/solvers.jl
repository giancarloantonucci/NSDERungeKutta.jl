"""
    IMEXEuler(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitExplicitRungeKuttaSolver
    IMEXSSP1_111(args...; kwargs...) :: ImplicitExplicitRungeKuttaSolver

returns an [`ImplicitExplicitRungeKuttaSolver`](@ref) for the 1st-order IMEXEuler method.
"""
function IMEXEuler(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    pᴵ = pᴱ = 1
    implicitableau = ButcherTableau(float([
         1 1;
        pᴵ 1;
    ]))
    explicitableau = ButcherTableau(float([
         0 0;
        pᴱ 1;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return IERK(implicitableau, explicitableau, h, newton)
end
@doc (@doc IMEXEuler) IMEXSSP1_111(args...; kwargs...) = IMEXEuler(args...; kwargs...)

"""
    IMEXSSP2_222(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitExplicitRungeKuttaSolver

returns an [`ImplicitExplicitRungeKuttaSolver`](@ref) for the 2nd-order IMEX-SSP2(2,2,2) L-stable scheme.
"""
function IMEXSSP2_222(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    pᴵ = pᴱ = 2
    γ = 1 - 1/√2
    implicitableau = ButcherTableau(float([
          γ    γ   0;
        1-γ 1-2γ   γ;
         pᴵ  1/2 1/2;
        ]))
    explicitableau = ButcherTableau(float([
         0   0   0;
         1   1   0;
        pᴱ 1/2 1/2;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return IERK(implicitableau, explicitableau, h, newton)
end

"""
    IMEXSSP2_322(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitExplicitRungeKuttaSolver

returns an [`ImplicitExplicitRungeKuttaSolver`](@ref) for the 2nd-order IMEX-SSP2(3,2,2) stiffly-accurate scheme.
"""
function IMEXSSP2_322(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    pᴵ = pᴱ = 2
    implicitableau = ButcherTableau(float([
        1/2  1/2   0   0;
          0 -1/2 1/2   0;
          1    0 1/2 1/2;
         pᴵ    0 1/2 1/2;
        ]))
    explicitableau = ButcherTableau(float([
         0   0   0   0;
         1   0   0   0;
         1   0   1   0;
        pᴱ   0 1/2 1/2;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return IERK(implicitableau, explicitableau, h, newton)
end

"""
    IMEXSSP2_332(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitExplicitRungeKuttaSolver

returns an [`ImplicitExplicitRungeKuttaSolver`](@ref) for the 2nd-order IMEX-SSP2(3,3,2) stiffly-accurate scheme.
"""
function IMEXSSP2_332(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    pᴵ = pᴱ = 2
    implicitableau = ButcherTableau(float([
        1/4 1/4   0   0;
        1/4   0 1/4   0;
          1 1/3 1/3 1/3;
         pᴵ 1/3 1/3 1/2;
    ]))
    explicitableau = ButcherTableau(float([
          0   0   0   0;
        1/2 1/2   0   0;
          1 1/2 1/2   0;
         pᴱ 1/3 1/3 1/3;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return IERK(implicitableau, explicitableau, h, newton)
end

"""
    IMEXSSP3_332(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10) :: ImplicitExplicitRungeKuttaSolver

returns an [`ImplicitExplicitRungeKuttaSolver`](@ref) for the 3rd-order IMEX-SSP3(3,3,2) L-stable scheme.
"""
function IMEXSSP3_332(; h::Real=0.0, rtol::Real=1e-3, nits::Integer=10)
    pᴵ = pᴱ = 3
    γ = 1 - 1/√2
    implicitableau = ButcherTableau(float([
          γ     γ   0   0;
        1-γ  1-2γ   γ   0;
        1/2 1/2-γ   0   γ;
         pᴵ   1/6 1/6 2/3;
    ]))
    explicitableau = ButcherTableau(float([
          0   0   0   0;
          1   1   0   0;
        1/2 1/4 1/4   0;
         pᴱ 1/6 1/6 2/3;
    ]))
    newton = NewtonParameters(rtol=rtol, nits=nits)
    return IERK(implicitableau, explicitableau, h, newton)
end
