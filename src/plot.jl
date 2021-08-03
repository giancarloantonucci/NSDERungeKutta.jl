@recipe function f(solution::RungeKuttaSolution; vars = nothing)
    fontfamily     --> "Computer Modern"
    framestyle     --> :box
    gridalpha      --> 0.2
    legend         --> :none
    linewidth      --> 1.5
    minorgrid      --> 0.1
    minorgridstyle --> :dash
    seriestype     --> :path
    xwiden         --> false
    tick_direction --> :out
    (L, N) = size(solution)
    if vars isa Tuple && length(vars) ≥ 1
        vars = vars
    elseif vars isa Integer
        vars = (vars, )
    elseif vars isa Nothing
        vars = 1:L
    else
        error("Got $(typeof(vars)) instead of `Tuple` with `length(vars)` ≥ 1 or `Integer`.")
    end
    @↓ u, t = solution
    [(t, [u[n][i] for n = 1:N]) for i in vars]
end

@userplot PHASEPLOT
@recipe function f(h::PHASEPLOT; vars = nothing)
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    solution = h.args[1] isa RungeKuttaSolution ? h.args[1] : error("Got $(typeof(h.args)) instead of `RungeKuttaSolution`.")
    fontfamily     --> "Computer Modern"
    framestyle     --> :box
    gridalpha      --> 0.2
    legend         --> :none
    linewidth      --> 1.5
    minorgrid      --> 0.1
    minorgridstyle --> :dash
    seriestype     --> :path
    tick_direction --> :out
    (L, N) = size(solution)
    if vars isa Tuple && length(vars) ≥ 1
        vars = vars
    elseif vars isa Integer
        vars = (vars, )
    elseif vars isa Nothing
        vars = 1:L
    else
        error("Got $(typeof(vars)) instead of `Tuple` with `length(vars)` ≥ 1 or `Integer`.")
    end
    @↓ u, t = solution
    tuple([[u[n][i] for n = 1:N] for i in vars]...)
end

@userplot STABILITY
@recipe function f(h::STABILITY; lims = (-4, 4), xlims = lims, ylims = lims)
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    R = h.args[1] isa ButcherTableau   ? z -> ℛ(z, h.args[1])         :
        h.args[1] isa RungeKuttaSolver ? z -> ℛ(z, h.args[1].tableau) :
        h.args[1] isa Function         ? h.args[1]                    :
        error("Got $(typeof(h.args)) instead of `ButcherTableau`, `RungeKuttaSolver` or `Function`.")
    fontfamily     --> "Computer Modern"
    framestyle     --> :box
    gridalpha      --> 0.2
    gridstyle      --> :dot
    legend         --> :none
    levels         --> [1.0]
    linewidth      --> 1.5
    minorgrid      --> 0.1
    minorgridstyle --> :dash
    seriestype     --> :contour
    tick_direction --> :out
    Δx = max(1, abs(xlims[2] - xlims[1]))
    Δy = max(1, abs(ylims[2] - ylims[1]))
    x = length(xlims) == 2 ? LinRange(xlims..., Int(Δx * 101)) : LinRange(xlims...)
    y = length(ylims) == 2 ? LinRange(ylims..., Int(Δy * 101)) : LinRange(ylims...)
    function f(x, y)
        z = x + 1im * y
        p = abs(R(z))
    end
    x, y, f
end

@userplot STABILITYF
@recipe function f(h::STABILITYF; lims = (-4, 4), xlims = lims, ylims = lims)
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    R = h.args[1] isa ButcherTableau   ? z -> ℛ(z, h.args[1])         :
        h.args[1] isa RungeKuttaSolver ? z -> ℛ(z, h.args[1].tableau) :
        h.args[1] isa Function         ? h.args[1]                    :
        error("Got $(typeof(h.args)) instead of `ButcherTableau`, `RungeKuttaSolver` or `Function`.")
    fontfamily     --> "Computer Modern"
    clims          --> (0, 1)
    colorbar       --> true
    framestyle     --> :box
    gridalpha      --> 0.2
    legend         --> :none
    minorgrid      --> 0.1
    minorgridstyle --> :dash
    seriestype     --> :heatmap
    tick_direction --> :out
    Δx = max(1, abs(xlims[2] - xlims[1]))
    Δy = max(1, abs(ylims[2] - ylims[1]))
    x = length(xlims) == 2 ? LinRange(xlims..., Int(Δx * 101)) : LinRange(xlims...)
    y = length(ylims) == 2 ? LinRange(ylims..., Int(Δy * 101)) : LinRange(ylims...)
    function f(x, y)
        z = x + 1im * y
        p = abs(R(z))
        (p > 1 ? -Inf : p)
    end
    x, y, f
end

@userplot ORDERSTAR
@recipe function f(h::ORDERSTAR; lims = (-4, 4), xlims = lims, ylims = lims)
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    R = h.args[1] isa ButcherTableau   ? z -> ℛ(z, h.args[1])         :
        h.args[1] isa RungeKuttaSolver ? z -> ℛ(z, h.args[1].tableau) :
        h.args[1] isa Function         ? h.args[1]                    :
        error("Got $(typeof(h.args)) instead of `ButcherTableau`, `RungeKuttaSolver` or `Function`.")
    fontfamily     --> "Computer Modern"
    framestyle     --> :box
    gridalpha      --> 0.2
    legend         --> :none
    levels         --> [1.0]
    linewidth      --> 1.5
    minorgrid      --> 0.1
    minorgridstyle --> :dash
    seriestype     --> :contour
    tick_direction --> :out
    Δx = max(1, abs(xlims[2] - xlims[1]))
    Δy = max(1, abs(ylims[2] - ylims[1]))
    x = length(xlims) == 2 ? LinRange(xlims..., Int(Δx * 101)) : LinRange(xlims...)
    y = length(ylims) == 2 ? LinRange(ylims..., Int(Δy * 101)) : LinRange(ylims...)
    function f(x, y)
        z = x + 1im * y
        abs(R(z) * exp(-z))
    end
    x, y, f
end

@userplot ORDERSTARF
@recipe function f(h::ORDERSTARF; lims = (-4, 4), xlims = lims, ylims = lims)
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    R = h.args[1] isa ButcherTableau   ? z -> ℛ(z, h.args[1])         :
        h.args[1] isa RungeKuttaSolver ? z -> ℛ(z, h.args[1].tableau) :
        h.args[1] isa Function         ? h.args[1]                    :
        error("Got $(typeof(h.args)) instead of `ButcherTableau`, `RungeKuttaSolver` or `Function`.")
    fontfamily     --> "Computer Modern"
    clims          --> (0, 1)
    colorbar       --> true
    framestyle     --> :box
    gridalpha      --> 0.2
    legend         --> :none
    minorgrid      --> 0.1
    minorgridstyle --> :dash
    seriestype     --> :heatmap
    tick_direction --> :out
    Δx = max(1, abs(xlims[2] - xlims[1]))
    Δy = max(1, abs(ylims[2] - ylims[1]))
    x = length(xlims) == 2 ? LinRange(xlims..., Int(Δx * 101)) : LinRange(xlims...)
    y = length(ylims) == 2 ? LinRange(ylims..., Int(Δy * 101)) : LinRange(ylims...)
    function f(x, y)
        z = x + 1im * y
        p = abs(R(z) * exp(-z))
        (p > 1 ? NaN : p)
    end
    x, y, f
end
