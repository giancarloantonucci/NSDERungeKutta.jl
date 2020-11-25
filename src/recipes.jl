@recipe function f(solution::RungeKuttaSolution; vars=nothing)
    framestyle --> :box
    legend     --> :none
    linewidth  --> 1.5
    seriestype --> :path
    @↓ u, t = solution
    L = length(u[Block(1)])
    N = length(t)
    if vars != nothing && !(vars isa Tuple && length(vars) > 1)
        error("Got $(typeof(vars)) instead of `Tuple` with length > 1.")
    elseif vars == nothing
        vars = 1:L
    end
    [(t, [u[BlockIndex(n, i)] for n = 1:N]) for i in vars]
end

@userplot PHASEPLOT
@recipe function f(h::PHASEPLOT; vars=nothing)
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    solution = h.args[1] isa RungeKuttaSolution ? h.args[1] : error("Got $(typeof(h.args)) instead of `RungeKuttaSolution`.")
    framestyle --> :box
    legend     --> :none
    seriestype --> :path
    @↓ u, t = solution
    L = length(u[Block(1)])
    N = length(t)
    if vars != nothing && !(vars isa Tuple && length(vars) > 1)
        error("Got $(typeof(vars)) instead of `Tuple` with length > 1.")
    elseif vars == nothing
        vars = 1:L
    end
    tuple([[u[BlockIndex(n, i)] for n = 1:N] for i in vars]...)
end

@userplot STABILITY
@recipe function f(h::STABILITY; xspan=(-3, 3), yspan=(-3, 3))
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    R = h.args[1] isa ButcherTableau   ? z -> ℛ(z, h.args[1])         :
        h.args[1] isa RungeKuttaSolver ? z -> ℛ(z, h.args[1].tableau) :
        h.args[1] isa Function         ? h.args[1]                    :
        error("Got $(typeof(h.args)) instead of `ButcherTableau`, `RungeKuttaSolver` or `Function`.")
    framestyle --> :box
    legend     --> :none
    levels     --> [1.0]
    seriestype --> :contour
    x = length(xspan) == 2 ? LinRange(xspan..., 501) : LinRange(xspan...)
    y = length(yspan) == 2 ? LinRange(yspan..., 501) : LinRange(yspan...)
    f(x, y) = abs(R(x + 1im*y))
    x, y, f
end

@userplot STABILITYF
@recipe function f(h::STABILITYF; xspan=(-3, 3), yspan=(-3, 3))
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    R = h.args[1] isa ButcherTableau   ? z -> ℛ(z, h.args[1])         :
        h.args[1] isa RungeKuttaSolver ? z -> ℛ(z, h.args[1].tableau) :
        h.args[1] isa Function         ? h.args[1]                    :
        error("Got $(typeof(h.args)) instead of `ButcherTableau`, `RungeKuttaSolver` or `Function`.")
    colorbar   --> :true
    framestyle --> :box
    legend     --> :none
    seriestype --> :heatmap
    x = length(xspan) == 2 ? LinRange(xspan..., 501) : LinRange(xspan...)
    y = length(yspan) == 2 ? LinRange(yspan..., 501) : LinRange(yspan...)
    function f(x, y)
        p = abs(R(x + 1im*y))
        (p > 1 ? NaN : p)
    end
    x, y, f
end

@userplot ORDERSTAR
@recipe function f(h::ORDERSTAR; xspan=(-3, 3), yspan=(-3, 3))
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    R = h.args[1] isa ButcherTableau   ? z -> ℛ(z, h.args[1])         :
        h.args[1] isa RungeKuttaSolver ? z -> ℛ(z, h.args[1].tableau) :
        h.args[1] isa Function         ? h.args[1]                    :
        error("Got $(typeof(h.args)) instead of `ButcherTableau`, `RungeKuttaSolver` or `Function`.")
    framestyle --> :box
    legend     --> :none
    levels     --> [1.0]
    seriestype --> :contour
    x = length(xspan) == 2 ? LinRange(xspan..., 501) : LinRange(xspan...)
    y = length(yspan) == 2 ? LinRange(yspan..., 501) : LinRange(yspan...)
    function f(x, y)
        z = x + 1im * y
        abs(R(z) * exp(-z))
    end
    x, y, f
end

@userplot ORDERSTARF
@recipe function f(h::ORDERSTARF; xspan=(-3, 3), yspan=(-3, 3))
    length(h.args) == 1 ? true : error("Got too many arguments: $(length(h.args)).")
    R = h.args[1] isa ButcherTableau   ? z -> ℛ(z, h.args[1])         :
        h.args[1] isa RungeKuttaSolver ? z -> ℛ(z, h.args[1].tableau) :
        h.args[1] isa Function         ? h.args[1]                    :
        error("Got $(typeof(h.args)) instead of `ButcherTableau`, `RungeKuttaSolver` or `Function`.")
    colorbar   --> :true
    framestyle --> :box
    legend     --> :none
    seriestype --> :heatmap
    x = length(xspan) == 2 ? LinRange(xspan..., 501) : LinRange(xspan...)
    y = length(yspan) == 2 ? LinRange(yspan..., 501) : LinRange(yspan...)
    function f(x, y)
        z = x + 1im * y
        p = abs(R(z) * exp(-z))
        (p > 1 ? NaN : p)
    end
    x, y, f
end
