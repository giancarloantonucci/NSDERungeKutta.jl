# RungeKutta.jl

A Julia package implementing Runge-Kutta methods.

[![Docs Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://giancarloantonucci.github.io/RungeKutta.jl/stable) [![Docs Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://giancarloantonucci.github.io/RungeKutta.jl/dev) [![Build Status](https://img.shields.io/github/workflow/status/giancarloantonucci/RungeKutta.jl/CI)](https://github.com/giancarloantonucci/RungeKutta.jl/actions) [![Coverage](https://img.shields.io/codecov/c/github/giancarloantonucci/RungeKutta.jl?label=coverage)](https://codecov.io/gh/giancarloantonucci/RungeKutta.jl)

## Installation

RungeKutta.jl is compatible with Julia v1.0 and above. From the Julia REPL,

```julia
]add https://github.com/giancarloantonucci/RungeKutta.jl
```

Below is a brief introduction, but read the [documentation](https://giancarloantonucci.github.io/RungeKutta.jl/dev) for a complete overview of this package.

## Usage

Let's say that we want to solve the [simple gravity pendulum problem](https://en.wikipedia.org/wiki/Pendulum_(mathematics)#Simple_gravity_pendulum) using the [midpoint method](https://en.wikipedia.org/wiki/Midpoint_method). Here is one way to do it with RungeKutta.jl:

```julia
using RungeKutta
f(u, t) = [u[2]; -9.81 * sin(u[1])]
u0 = [0.0; π/2]
tspan = (0.0, 2π)
problem = IVP(f, u0, tspan)
solver = Midpoint(h = 1e-2)
solution = solve(problem, solver)
```

We can plot the obtained `solution` by extracting its fields `u` and `t`, e.g. using the convenient macro `@↓ u, t = solution` from [ArrowMacros.jl](https://github.com/giancarloantonucci/ArrowMacros.jl). Alternatively, we can use the available predefined recipes:

```julia
using Plots, LaTeXStrings
p₁ = plot(solution, xlabel = L"t", label = [L"\theta" L"\omega"], legend = true)
p₂ = phaseplot(solution, vars = (1, 2), xlabel = L"\theta", ylabel = L"\omega")
plot(size = (900, 450), p₁, p₂)
```

![svg](imgs/pendulum.svg)

For convenience, RungeKutta.jl re-exports all the ODE problems defined in [NSDEBase.jl](https://github.com/giancarloantonucci/NSDEBase.jl), e.g. `SimplePendulum` for the above problem.

RungeKutta.jl has some predefined recipes to plot **stability regions** and **order stars** too:

```julia
p₁ = stabilityf(RK4(), xlabel = L"\Re(z)", ylabel = L"\Im(z)", colour = :blues)
p₂ = orderstarf(RK4(), xlabel = L"\Re(z)", ylabel = L"\Im(z)", colour = :blues)
plot(size = (1000, 400), p₁, p₂, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
```

![svg](imgs/regions.svg)

## Available methods

RungeKutta.jl currently supports the following methods:

**Explicit**: `Euler`/`ExplicitEuler`, `Midpoint`/`ExplicitMidpoint`, `Heun2`, `Ralston2`, `Heun3`, `Kutta3`, `Ralston3`, `SSPRK3`, `RK4`, `Rule38`, `HeunEuler`, `Fehlberg45`/`F45`, `DormandPrince54`/`DP54`, `Verner65`/`V65`.

**Implicit**: `BackwardEuler`/`ImplicitEuler`, `ImplicitMidpoint`, `CrankNicolson`, `SDIRK3`, `GaussLegendre4`/`GL4`, `GaussLegendre6`/`GL6`, `LobattoIIIA4`, `LobattoIIIB2`, `LobattoIIIB4`, `LobattoIIIC2`, `LobattoIIIC4`, `RadauIA3`, `RadauIA5`, `RadauIIA3`, `RadauIIA5`.
