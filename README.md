# RungeKutta

A Julia package implementing Runge-Kutta methods.

[![Docs Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://antonuccig.github.io/RungeKutta.jl/stable) [![Docs Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://antonuccig.github.io/RungeKutta.jl/dev) [![Build Status](https://img.shields.io/github/workflow/status/antonuccig/RungeKutta.jl/CI)](https://github.com/antonuccig/RungeKutta.jl/actions) [![Coverage](https://img.shields.io/codecov/c/github/antonuccig/RungeKutta.jl?label=codecov)](https://codecov.io/gh/antonuccig/RungeKutta.jl)

## Installation

`RungeKutta` is compatible with Julia `v1.0` and above. From the Julia REPL,

```julia
]add https://github.com/antonuccig/RungeKutta.jl
```

## Usage

Let's say that we want to solve the [simple gravity pendulum problem](https://en.wikipedia.org/wiki/Pendulum_(mathematics)#Simple_gravity_pendulum) using the [midpoint method](https://en.wikipedia.org/wiki/Midpoint_method). Here is how to do it with `RungeKutta`:

```julia
using RungeKutta
f(u, t) = [u[2]; -9.81 * sin(u[1])]
u0 = [0.0; π/2]
tspan = (0.0, 2π)
problem = IVP(f, u0, tspan)
solver = Midpoint(h = 1e-2)
solution = solve(problem, solver)
```

### Plotting

We can plot the obtained `solution` by extracting its fields `u` and `t`, e.g. with the convenient macro `@↓ u, t = solution` from `ArrowMacros`. Alternatively, we can use the predefined recipes:

```julia
using Plots, LaTeXStrings
default(fontfamily = "Computer Modern")
p₁ = plot(solution, xlabel = L"t", label = [L"\theta" L"\omega"], legend = true)
p₂ = phaseplot(solution, vars = (1, 2), xlabel = L"\theta", ylabel = L"\omega")
plot(size = (800, 400), p₁, p₂)
# savefig("pendulum.svg")
```

![svg](images/pendulum.svg)

### Predefined ODE recipes

`RungeKutta` comes with some predefined ODE problems, like the [Lorenz system](https://en.wikipedia.org/wiki/Lorenz_system):

```julia
u0 = [2.0, 3.0, -14.0]
tspan = (0.0, 10.0)
problem = Lorenz(u0, tspan)
solver = F45(h = 1e-3)
solution = solve(problem, solver)
plot(solution, xlabel = L"t", label = [L"x" L"y" L"z"], legend = true)
# savefig("lorenz.svg")
```

![svg](images/lorenz.svg)

### Stability regions and order stars

`RungeKutta` has also predefined recipes to plot stability regions and order stars:

```julia
p₁ = stabilityf(RK4(), xlabel = L"\Re(z)", ylabel = L"\Im(z)", colour = :blues)
p₂ = orderstarf(RK4(), xlabel = L"\Re(z)", ylabel = L"\Im(z)", colour = :blues)
plot(size = (1000, 400), p₁, p₂, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
# savefig("regions.svg")
```

![svg](images/regions.svg)

## Available methods

`RungeKutta` currently supports the following methods:

<!-- explicit (`Euler`/`ExplicitEuler`, `Midpoint`/`ExplicitMidpoint`, `Heun2`, `Ralston2`, `Heun3`, `Kutta3`, `Ralston3`, `SSPRK3`, `RK4`, `Rule38`, `HeunEuler`, `Fehlberg45`/`F45`, `DormandPrince54`/`DP54`, `Verner65`/`V65`) and implicit methods (`BackwardEuler`/`ImplicitEuler`, `ImplicitMidpoint`, `CrankNicolson`, `SDIRK3`, `GaussLegendre4`/`GL4`, `GaussLegendre6`/`GL6`, `LobattoIIIA4`, `LobattoIIIB2`, `LobattoIIIB4`, `LobattoIIIC2`, `LobattoIIIC4`, `RadauIA3`, `RadauIA5`, `RadauIIA3`, `RadauIIA5`). -->

<details><summary>Explicit</summary>

- <code>Euler</code>/<code>ExplicitEuler</code>
- <code>Midpoint</code>/<code>ExplicitMidpoint</code>
- <code>Heun2</code>
- <code>Ralston2</code>
- <code>Heun3</code>
- <code>Kutta3</code>
- <code>Ralston3</code>
- <code>SSPRK3</code>
- <code>RK4</code>
- <code>Rule38</code>
- <code>HeunEuler</code>
- <code>Fehlberg45</code>/<code>F45</code>
- <code>DormandPrince54</code>/<code>DP54</code>
- <code>Verner65</code>/<code>V65</code></details>



<details><summary>Implicit</summary>

- <code>BackwardEuler</code>/<code>ImplicitEuler</code>
- <code>ImplicitMidpoint</code>
- <code>CrankNicolson</code>
- <code>SDIRK3</code>
- <code>GaussLegendre4</code>/<code>GL4</code>
- <code>GaussLegendre6</code>/<code>GL6</code>
- <code>LobattoIIIA4</code>
- <code>LobattoIIIB2</code>
- <code>LobattoIIIB4</code>
- <code>LobattoIIIC2</code>
- <code>LobattoIIIC4</code>
- <code>RadauIA3</code>
- <code>RadauIA5</code>
- <code>RadauIIA3</code>
- <code>RadauIIA5</code></details>

## What's next?

Current plans for future developments are:

- Improve performance and error messages.
- Automatic size detection of stability region.
- IMEX methods.
