@recipe function f(solution::RungeKuttaSolution; variables=nothing, is_complex=false)
    gridalpha --> 0.2
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :path
    (d, N) = size(solution)
    if variables isa Nothing
        variables = 1:d
    elseif variables isa Integer
        variables = [variables]
    elseif variables isa Tuple
        # variables = [variable for variable in variables]
    elseif variables isa AbstractVector
        # variables = variables
    else
        throw(ArgumentError("Got $(typeof(variables)) instead of either `Integer`, `Tuple`, or `AbstractVector`."))
    end
    @↓ u, t = solution
    if is_complex
        for i in variables
            @series begin
                seriescolor --> i
                ([real.(u[n][i]) for n = 1:N], [imag.(u[n][i]) for n = 1:N])
                # (t, [real.(u[n][i]) for n = 1:N], [imag.(u[n][i]) for n = 1:N])
            end
        end
    else
        for (i, variable) in enumerate(variables)
            @series begin
                xwiden --> false
                if haskey(plotattributes, :label) && plotattributes[:label] isa AbstractVector
                    label := plotattributes[:label][i]
                end
                seriescolor --> i
                (t, [u[n][variable] for n = 1:N])
            end
        end
    end
end

@recipe function f(wrapper::NSDEBase._PhasePlot{<:RungeKuttaSolution})
    gridalpha --> 0.2
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :path
    @↓ solution ← plottable = wrapper
    (d, N) = size(solution)
    variables = haskey(plotattributes, :variables) ? plotattributes[:variables] : (1:d)
    @↓ u, t = solution
    return tuple([[u[n][i] for n = 1:N] for i in variables]...)
end

@userplot STABILITY
@recipe function f(h::STABILITY; resolution=100, span=range(-5,5,length=resolution), xspan=span, yspan=span)
    R = if h.args[1] isa ButcherTableau
        z -> ℛ(z, h.args[1])
    elseif h.args[1] isa AbstractRungeKuttaSolver
        z -> ℛ(z, h.args[1].tableau)
    elseif h.args[1] isa Function
        h.args[1]
    else
        throw(ArgumentError("Got $(typeof(h.args)) instead of either `ButcherTableau`, `AbstractRungeKuttaSolver`, or `Function`."))
    end
    legend --> false
    levels --> [1.0]
    seriestype --> :contour
    function f(x, y)
        z = x + 1im * y
        return abs(R(z))
    end
    return xspan, yspan, f
end

@userplot STABILITYF
@recipe function f(h::STABILITYF; resolution=100, span=range(-5,5,length=resolution), xspan=span, yspan=span)
    R = if h.args[1] isa ButcherTableau
        z -> ℛ(z, h.args[1])
    elseif h.args[1] isa AbstractRungeKuttaSolver
        z -> ℛ(z, h.args[1].tableau)
    elseif h.args[1] isa Function
        h.args[1]
    else
        throw(ArgumentError("Got $(typeof(h.args)) instead of either `ButcherTableau`, `AbstractRungeKuttaSolver`, or `Function`."))
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
@recipe function f(h::ORDERSTAR; resolution=100, span=range(-5,5,length=resolution), xspan=span, yspan=span)
    R = if h.args[1] isa ButcherTableau
        z -> ℛ(z, h.args[1])
    elseif h.args[1] isa AbstractRungeKuttaSolver
        z -> ℛ(z, h.args[1].tableau)
    elseif h.args[1] isa Function
        h.args[1]
    else
        throw(ArgumentError("Got $(typeof(h.args)) instead of either `ButcherTableau`, `AbstractRungeKuttaSolver`, or `Function`."))
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
@recipe function f(h::ORDERSTARF; resolution=100, span=range(-5,5,length=resolution), xspan=span, yspan=span)
    R = if h.args[1] isa ButcherTableau
        z -> ℛ(z, h.args[1])
    elseif h.args[1] isa AbstractRungeKuttaSolver
        z -> ℛ(z, h.args[1].tableau)
    elseif h.args[1] isa Function
        h.args[1]
    else
        throw(ArgumentError("Got $(typeof(h.args)) instead of either `ButcherTableau`, `AbstractRungeKuttaSolver`, or `Function`."))
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

@recipe function f(ts::AbstractVector, hs::StepSizes)
    gridalpha --> 0.2
    @series begin
        markershape --> :circle
        markerstrokewidth --> 0
        return ts, hs.accepted
    end
    for n in 1:length(hs.rejected)
        @series begin
            markershape --> :x
            seriestype --> :scatter
            seriescolor --> 2
            y = hs.rejected[n]
            x = ts[n] * ones(length(y))
            return x, y
        end
    end
end

@recipe function f(hs::StepSizes)
    gridalpha --> 0.2
    @series begin
        markershape --> :circle
        markerstrokewidth --> 0
        return hs.accepted
    end
    for n in 1:length(hs.rejected)
        @series begin
            markershape --> :x
            seriestype --> :scatter
            seriescolor --> 2
            y = hs.rejected[n]
            x = n * ones(length(y))
            return x, y
        end
    end
end
