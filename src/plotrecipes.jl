@recipe function f(solution::RungeKuttaSolution; vars=nothing, iscomplex=false, plotkind=:plot, label="")
    framestyle --> :box
    gridalpha --> 0.2
    linewidth --> 1.5
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :path
    tick_direction --> :out
    (L, N) = size(solution)
    if vars isa AbstractVector && length(vars) ≥ 1
        vars = vars
    elseif vars isa Integer
        vars = [vars]
    elseif vars isa Nothing
        vars = 1:L
    else
        error("Got $(typeof(vars)) instead of `AbstractVector` with `length(vars)` ≥ 1 or `Integer`.")
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
        if plotkind == :phaseplot
            @series begin
                seriescolor --> 1
                tuple([[u[n][i] for n = 1:N] for i in vars]...)
            end
        else
            xwiden --> false
            for (i, var) in enumerate(vars)
                @series begin
                    label := label isa AbstractVector ? label[i] : label
                    seriescolor --> i
                    (t, [u[n][var] for n = 1:N])
                end
            end
        end
    end
end

@userplot STABILITY
@recipe function f(h::STABILITY; lims=(-4, 4), xlims=lims, ylims=lims)
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    R = h.args[1] isa ButcherTableau ? z -> ℛ(z, h.args[1]) :
        h.args[1] isa AbstractRungeKuttaSolver ? z -> ℛ(z, h.args[1].tableau) :
        h.args[1] isa Function ? h.args[1] :
        error("Got $(typeof(h.args)) instead of `ButcherTableau`, `AbstractRungeKuttaSolver` or `Function`.")
    framestyle --> :box
    gridalpha --> 0.2
    gridstyle --> :dot
    legend --> false
    levels --> [1.0]
    linewidth --> 1.5
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :contour
    tick_direction --> :out
    Δx = max(1, abs(xlims[2] - xlims[1]))
    Δy = max(1, abs(ylims[2] - ylims[1]))
    x = length(xlims) == 2 ? LinRange(xlims..., Int(Δx * 101)) : LinRange(xlims...)
    y = length(ylims) == 2 ? LinRange(ylims..., Int(Δy * 101)) : LinRange(ylims...)
    function f(x, y)
        z = x + 1im * y
        return p = abs(R(z))
    end
    return x, y, f
end

@userplot STABILITYF
@recipe function f(h::STABILITYF; lims=(-4, 4), xlims=lims, ylims=lims)
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    R = h.args[1] isa ButcherTableau ? z -> ℛ(z, h.args[1]) :
        h.args[1] isa AbstractRungeKuttaSolver ? z -> ℛ(z, h.args[1].tableau) :
        h.args[1] isa Function ? h.args[1] :
        error("Got $(typeof(h.args)) instead of `ButcherTableau`, `AbstractRungeKuttaSolver` or `Function`.")
    clims --> (0, 1)
    colorbar --> true
    framestyle --> :box
    gridalpha --> 0.2
    legend --> false
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :heatmap
    tick_direction --> :out
    Δx = max(1, abs(xlims[2] - xlims[1]))
    Δy = max(1, abs(ylims[2] - ylims[1]))
    x = length(xlims) == 2 ? LinRange(xlims..., Int(Δx * 101)) : LinRange(xlims...)
    y = length(ylims) == 2 ? LinRange(ylims..., Int(Δy * 101)) : LinRange(ylims...)
    function f(x, y)
        z = x + 1im * y
        p = abs(R(z))
        return (p > 1 ? -Inf : p)
    end
    return x, y, f
end

@userplot ORDERSTAR
@recipe function f(h::ORDERSTAR; lims=(-4, 4), xlims=lims, ylims=lims)
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    R = h.args[1] isa ButcherTableau ? z -> ℛ(z, h.args[1]) :
        h.args[1] isa AbstractRungeKuttaSolver ? z -> ℛ(z, h.args[1].tableau) :
        h.args[1] isa Function ? h.args[1] :
        error("Got $(typeof(h.args)) instead of `ButcherTableau`, `AbstractRungeKuttaSolver` or `Function`.")
    framestyle --> :box
    gridalpha --> 0.2
    legend --> false
    levels --> [1.0]
    linewidth --> 1.5
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :contour
    tick_direction --> :out
    Δx = max(1, abs(xlims[2] - xlims[1]))
    Δy = max(1, abs(ylims[2] - ylims[1]))
    x = length(xlims) == 2 ? LinRange(xlims..., Int(Δx * 101)) : LinRange(xlims...)
    y = length(ylims) == 2 ? LinRange(ylims..., Int(Δy * 101)) : LinRange(ylims...)
    function f(x, y)
        z = x + 1im * y
        return abs(R(z) * exp(-z))
    end
    return x, y, f
end

@userplot ORDERSTARF
@recipe function f(h::ORDERSTARF; lims=(-4, 4), xlims=lims, ylims=lims)
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    R = h.args[1] isa ButcherTableau ? z -> ℛ(z, h.args[1]) :
        h.args[1] isa AbstractRungeKuttaSolver ? z -> ℛ(z, h.args[1].tableau) :
        h.args[1] isa Function ? h.args[1] :
        error("Got $(typeof(h.args)) instead of `ButcherTableau`, `AbstractRungeKuttaSolver` or `Function`.")
    clims --> (0, 1)
    colorbar --> true
    framestyle --> :box
    gridalpha --> 0.2
    legend --> false
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :heatmap
    tick_direction --> :out
    Δx = max(1, abs(xlims[2] - xlims[1]))
    Δy = max(1, abs(ylims[2] - ylims[1]))
    x = length(xlims) == 2 ? LinRange(xlims..., Int(Δx * 101)) : LinRange(xlims...)
    y = length(ylims) == 2 ? LinRange(ylims..., Int(Δy * 101)) : LinRange(ylims...)
    function f(x, y)
        z = x + 1im * y
        p = abs(R(z) * exp(-z))
        return (p > 1 ? NaN : p)
    end
    return x, y, f
end
