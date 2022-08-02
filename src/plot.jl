@recipe function f(solution::RungeKuttaSolution; vars=nothing, iscomplex=false, plotkind=:plot)
    gridalpha --> 0.2
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :path
    (L, N) = size(solution)
    if vars isa Nothing
        vars = 1:L
    elseif vars isa Integer
        vars = [vars]
    elseif vars isa Tuple
        # vars = [var for var in vars]
    elseif vars isa AbstractVector
        # vars = vars
    else
        throw(ArgumentError("Got $(typeof(vars)) instead of `Union{Integer, Tuple, AbstractVector}` with `length(vars)` ≥ 1."))
    end
    @↓ u, t = solution
    if iscomplex
        for i in vars
            @series begin
                seriescolor --> i
                [u[n][i] for n = 1:N]
            end
        end
    else
        for (i, var) in enumerate(vars)
            @series begin
                if haskey(plotattributes, :label) && plotattributes[:label] isa AbstractVector
                    label := plotattributes[:label][i]
                end
                seriescolor --> i
                (t, [u[n][var] for n = 1:N])
            end
        end
    end
end

@recipe function f(wrappedobject::NSDEBase._PhasePlot{<:RungeKuttaSolution})
    gridalpha --> 0.2
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :path
    @↓ solution ← object = wrappedobject
    (L, N) = size(solution)
    vars = haskey(plotattributes, :vars) ? plotattributes[:vars] : (1:L)
    @↓ u, t = solution
    return tuple([[u[n][i] for n = 1:N] for i in vars]...)
end

@userplot STABILITY
@recipe function f(h::STABILITY; span=range(-5,5,length=100), xspan=span, yspan=span)
    R = if h.args[1] isa ButcherTableau
        z -> ℛ(z, h.args[1])
    elseif h.args[1] isa AbstractRungeKuttaSolver
        z -> ℛ(z, h.args[1].tableau)
    elseif h.args[1] isa Function
        h.args[1]
    else
        throw(ArgumentError("Got $(typeof(h.args)) instead of `Union{ButcherTableau, AbstractRungeKuttaSolver, Function}`."))
    end
    legend --> false
    levels --> [1.0]
    seriestype --> :contour
    function f(x, y)
        z = x + 1im * y
        return p = abs(R(z))
    end
    return xspan, yspan, f
end

@userplot STABILITYF
@recipe function f(h::STABILITYF; span=range(-5,5,length=100), xspan=span, yspan=span)
    R = if h.args[1] isa ButcherTableau
        z -> ℛ(z, h.args[1])
    elseif h.args[1] isa AbstractRungeKuttaSolver
        z -> ℛ(z, h.args[1].tableau)
    elseif h.args[1] isa Function
        h.args[1]
    else
        throw(ArgumentError("Got $(typeof(h.args)) instead of `Union{ButcherTableau, AbstractRungeKuttaSolver, Function}`."))
    end
    clims --> (0, 1)
    colorbar --> true
    legend --> false
    seriestype --> :heatmap
    function f(x, y)
        z = x + 1im * y
        p = abs(R(z))
        return (p > 1 ? -Inf : p)
    end
    return xspan, yspan, f
end

@userplot ORDERSTAR
@recipe function f(h::ORDERSTAR; span=range(-5,5,length=100), xspan=span, yspan=span)
    R = if h.args[1] isa ButcherTableau
        z -> ℛ(z, h.args[1])
    elseif h.args[1] isa AbstractRungeKuttaSolver
        z -> ℛ(z, h.args[1].tableau)
    elseif h.args[1] isa Function
        h.args[1]
    else
        throw(ArgumentError("Got $(typeof(h.args)) instead of `Union{ButcherTableau, AbstractRungeKuttaSolver, Function}`."))
    end
    legend --> false
    levels --> [1.0]
    seriestype --> :contour
    function f(x, y)
        z = x + 1im * y
        return abs(R(z) * exp(-z))
    end
    return xspan, yspan, f
end

@userplot ORDERSTARF
@recipe function f(h::ORDERSTARF; span=range(-5,5,length=100), xspan=span, yspan=span)
    R = if h.args[1] isa ButcherTableau
        z -> ℛ(z, h.args[1])
    elseif h.args[1] isa AbstractRungeKuttaSolver
        z -> ℛ(z, h.args[1].tableau)
    elseif h.args[1] isa Function
        h.args[1]
    else
        throw(ArgumentError("Got $(typeof(h.args)) instead of `Union{ButcherTableau, AbstractRungeKuttaSolver, Function}`."))
    end
    clims --> (0, 1)
    colorbar --> true
    legend --> false
    seriestype --> :heatmap
    function f(x, y)
        z = x + 1im * y
        p = abs(R(z) * exp(-z))
        return (p > 1 ? NaN : p)
    end
    return xspan, yspan, f
end
