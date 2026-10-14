stability_plot_input(input) = input
function stability_plot_input(input::AbstractRungeKuttaSolver)
    hasproperty(input, :tableau) && input.tableau isa ButcherTableau || throw(
        ArgumentError(
            "scalar stability plots support ordinary Runge-Kutta tableaus, not IMEX or exponential solvers",
        ),
    )
    return input.tableau
end

# The scalar function R(z) a stability or order-star recipe draws, from a method
# constructor, a configured solver, a tableau or R itself.
function stability_plot_function(input)
    input = stability_plot_input(input)
    input isa ButcherTableau && return z -> stability_function(z, input)
    input isa Function && return input
    return throw(
        ArgumentError(
            "Expected ButcherTableau, AbstractRungeKuttaSolver, or Function. Got $(typeof(input)) instead.",
        ),
    )
end

for method in (
    Euler,
    Heun2,
    Midpoint,
    Ralston2,
    Heun3,
    RungeKutta3,
    Ralston3,
    SSPRK3,
    Ralston4,
    RungeKutta4,
    Rule38,
    Butcher5,
    KuttaNystrom5,
    Butcher6,
    Butcher7,
    HeunEuler,
    BogackiShampine,
    Fehlberg45,
    DormandPrince54,
    Verner65,
    Fehlberg78,
    ExplicitEuler,
    ExplicitMidpoint,
    RK3,
    RK4,
    F45,
    DP54,
    V65,
    F78,
    BackwardEuler,
    ImplicitMidpoint,
    SDIRK2,
    LobattoIII2,
    CrankNicolson,
    SDIRK3,
    RadauI3,
    RadauII3,
    SDIRK4,
    LobattoIII4,
    ImplicitEuler,
    GaussLegendre2,
    LobattoIIIA2,
    LobattoIIIC2,
    RadauIA3,
    RadauIIA3,
    GaussLegendre4,
    LobattoIIIA4,
    LobattoIIIB4,
    LobattoIIIC4,
    RadauI5,
    RadauIA5,
    RadauII5,
    RadauIIA5,
    GaussLegendre6,
)
    @eval stability_plot_input(::typeof($method)) = $method(; h=1.0).tableau
end

RecipesBase.@recipe function f(
    solution::RungeKuttaSolution; variables=nothing, iscomplex=false, skip=1
)
    gridalpha --> 0.2
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :path

    (d, N) = size(solution)

    # Decide which components to draw
    idxs = if variables isa Nothing
        1:d
    elseif variables isa Integer
        [variables]
    elseif variables isa AbstractVector || variables isa Tuple
        variables
    else
        throw(
            ArgumentError(
                "Variables must be Integer, Tuple, or AbstractVector. Got $(typeof(variables))."
            ),
        )
    end

    @↓ u, t = solution

    if iscomplex
        # Trace paths in the complex plane
        for (j, i) in enumerate(idxs)
            RecipesBase.@series begin
                seriescolor --> i
                ([real.(u[n][i]) for n in 1:skip:N], [imag.(u[n][i]) for n in 1:skip:N])
            end
        end
    else
        # Draw component vs time
        for (j, i) in enumerate(idxs)
            RecipesBase.@series begin
                xwiden --> false

                # Override label if provided in plot command
                if haskey(plotattributes, :label) &&
                    plotattributes[:label] isa Union{AbstractArray,Tuple}
                    label := plotattributes[:label][j]
                end

                seriescolor --> i
                (t[1:skip:N], [u[n][i] for n in 1:skip:N])
            end
        end
    end
end

RecipesBase.@recipe function f(
    wrapper::NSDEBase._PhasePlot{<:RungeKuttaSolution}; skip=1, variables=nothing
)
    gridalpha --> 0.2
    minorgrid --> 0.1
    minorgridstyle --> :dash
    seriestype --> :path

    @↓ solution ← plottable = wrapper
    (d, N) = size(solution)

    # Pick default variables if none given
    vars = if variables === nothing
        (1:min(d, 3))
    elseif variables isa Integer
        (variables,)
    else
        variables
    end

    @↓ u, t = solution

    # Pack chosen components
    return Tuple([u[n][v] for n in 1:skip:N] for v in vars)
end

RecipesBase.@userplot STABILITY
RecipesBase.recipetype(::Val{:stability}, args...) = STABILITY(args)
RecipesBase.@recipe function f(
    h::STABILITY; resolution=100, span=range(-5, 5, length=resolution), xspan=span, yspan=span
)
    R = stability_plot_function(h.args[1])
    legend --> false
    levels --> [1.0]
    seriestype --> :contour
    function f(x, y)
        z = x + 1im * y
        return abs(R(z))
    end
    return xspan, yspan, f
end

RecipesBase.@userplot STABILITYF
RecipesBase.recipetype(::Val{:stabilityf}, args...) = STABILITYF(args)
RecipesBase.@recipe function f(
    h::STABILITYF; resolution=100, span=range(-5, 5, length=resolution), xspan=span, yspan=span
)
    R = stability_plot_function(h.args[1])
    clims --> (0, 1)
    colorbar --> true
    legend --> false
    seriestype --> :heatmap
    function f(x, y)
        z = x + 1im * y
        p = abs(R(z))
        return (p > 1 ? NaN : p)
    end
    return xspan, yspan, f
end

RecipesBase.@userplot ORDERSTAR
RecipesBase.recipetype(::Val{:orderstar}, args...) = ORDERSTAR(args)
RecipesBase.@recipe function f(
    h::ORDERSTAR; resolution=100, span=range(-5, 5, length=resolution), xspan=span, yspan=span
)
    R = stability_plot_function(h.args[1])
    legend --> false
    levels --> [1.0]
    seriestype --> :contour
    function f(x, y)
        z = x + 1im * y
        return abs(R(z) * exp(-z))
    end
    return xspan, yspan, f
end

RecipesBase.@userplot ORDERSTARF
RecipesBase.recipetype(::Val{:orderstarf}, args...) = ORDERSTARF(args)
RecipesBase.@recipe function f(
    h::ORDERSTARF; resolution=100, span=range(-5, 5, length=resolution), xspan=span, yspan=span
)
    R = stability_plot_function(h.args[1])
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

RecipesBase.@recipe function f(ts::AbstractVector, hs::StepSizes)
    gridalpha --> 0.2

    # Accepted steps
    RecipesBase.@series begin
        markershape --> :circle
        markerstrokewidth --> 0
        return ts, hs.accepted
    end

    # Rejected steps
    if !isempty(hs.rejected)
        x_rej = eltype(ts)[]
        y_rej = eltype(hs.accepted)[]

        for (n, rejections) in enumerate(hs.rejected)
            n > length(ts) && continue # the trailing (empty) bin has no abscissa
            for h in rejections
                push!(x_rej, ts[n])
                push!(y_rej, h)
            end
        end

        if !isempty(x_rej)
            RecipesBase.@series begin
                markershape --> :x
                seriestype --> :scatter
                seriescolor --> 2
                return x_rej, y_rej
            end
        end
    end
end

RecipesBase.@recipe function f(hs::StepSizes)
    gridalpha --> 0.2

    # Accepted steps
    RecipesBase.@series begin
        markershape --> :circle
        markerstrokewidth --> 0
        return hs.accepted
    end

    # Rejected steps
    if !isempty(hs.rejected)
        x_rej = Int[]
        y_rej = eltype(hs.accepted)[]

        for (n, rejections) in enumerate(hs.rejected)
            for h in rejections
                push!(x_rej, n)
                push!(y_rej, h)
            end
        end

        if !isempty(x_rej)
            RecipesBase.@series begin
                markershape --> :x
                seriestype --> :scatter
                seriescolor --> 2
                return x_rej, y_rej
            end
        end
    end
end
