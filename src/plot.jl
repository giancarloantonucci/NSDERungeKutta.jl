@recipe function f(solution::RungeKuttaSolution; variables=nothing, iscomplex=false, plotkind=:plot)
    gridalpha --> 0.2
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :path
    xwiden --> false
    (num_variables, num_time_steps) = size(solution)
    if variables isa Nothing
        variables = 1:num_variables
    elseif variables isa Integer
        variables = [variables]
    elseif variables isa Tuple
        # variables = [variable for variable in variables]
    elseif variables isa AbstractVector
        # variables = variables
    else
        throw(ArgumentError("Got $(typeof(variables)) instead of `Union{Integer, Tuple, AbstractVector}` with `length(variables)` ≥ 1."))
    end
    @↓ solution_variables ← u, solution_time_steps ← t = solution
    if iscomplex
        for i in variables
            @series begin
                seriescolor --> i
                [solution_variables[n][i] for n = 1:num_time_steps]
            end
        end
    else
        for (i, variable) in enumerate(variables)
            @series begin
                if haskey(plotattributes, :label) && plotattributes[:label] isa AbstractVector
                    label := plotattributes[:label][i]
                end
                seriescolor --> i
                (solution_time_steps, [solution_variables[n][variable] for n = 1:num_time_steps])
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
    (num_variables, num_time_steps) = size(solution)
    variables = haskey(plotattributes, :variables) ? plotattributes[:variables] : (1:num_variables)
    @↓ solution_variables ← u, solution_time_steps ← t = solution
    return tuple([[solution_variables[n][i] for n = 1:num_time_steps] for i in variables]...)
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
        throw(ArgumentError("Got $(typeof(h.args)) instead of `Union{ButcherTableau, AbstractRungeKuttaSolver, Function}`."))
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
@recipe function f(h::ORDERSTAR; resolution=100, span=range(-5,5,length=resolution), xspan=span, yspan=span)
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
@recipe function f(h::ORDERSTARF; resolution=100, span=range(-5,5,length=resolution), xspan=span, yspan=span)
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
