# RungeKutta.jl

A Julia package implementing Runge-Kutta methods.

[![Docs Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://antonuccig.github.io/RungeKutta.jl/stable) [![Docs Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://antonuccig.github.io/RungeKutta.jl/dev) [![Build Status](https://img.shields.io/github/workflow/status/antonuccig/RungeKutta.jl/CI)](https://github.com/antonuccig/RungeKutta.jl/actions) [![Coverage](https://img.shields.io/codecov/c/github/antonuccig/RungeKutta.jl?label=coverage)](https://codecov.io/gh/antonuccig/RungeKutta.jl)

## Installation

RungeKutta is compatible with Julia v1.0 and above. From the Julia REPL,

```julia
]add https://github.com/antonuccig/RungeKutta.jl
```

## Usage

Let's say that we want to solve the [simple gravity pendulum problem](https://en.wikipedia.org/wiki/Pendulum_(mathematics)#Simple_gravity_pendulum) using the [midpoint method](https://en.wikipedia.org/wiki/Midpoint_method). Here is how to do it with RungeKutta:

```julia
using RungeKutta
f(u, t) = [u[2]; -9.81 * sin(u[1])]
u0 = [0.0; π/2]
tspan = (0.0, 2π)
problem = IVP(f, u0, tspan)
solver = Midpoint(h = 1e-2)
solution = solve(problem, solver)
```

We can plot the obtained `solution` by extracting its fields `u` and `t`, e.g. with the convenient macro `@↓ u, t = solution` from [ArrowMacros.jl](https://github.com/antonuccig/ArrowMacros.jl). Alternatively, we can use the available predefined recipes:

```julia
using Plots, LaTeXStrings
default(fontfamily = "Computer Modern")
p₁ = plot(solution, xlabel = L"t", label = [L"\theta" L"\omega"], legend = true)
p₂ = phaseplot(solution, vars = (1, 2), xlabel = L"\theta", ylabel = L"\omega")
plot(size = (800, 400), p₁, p₂)
# savefig("pendulum.svg")
```

![svg](images/pendulum.svg)

For convenience, RungeKutta re-exports all ODE problems (pre-)defined in [NSDEBase.jl](https://github.com/antonuccig/NSDEBase.jl), e.g. `Lorenz` for [Lorenz's ODEs](https://en.wikipedia.org/wiki/Lorenz_system):

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

RungeKutta also has some predefined recipes to plot stability regions and order stars:

```julia
p₁ = stabilityf(RK4(), xlabel = L"\Re(z)", ylabel = L"\Im(z)", colour = :blues)
p₂ = orderstarf(RK4(), xlabel = L"\Re(z)", ylabel = L"\Im(z)", colour = :blues)
plot(size = (1000, 400), p₁, p₂, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
# savefig("regions.svg")
```

![svg](images/regions.svg)

## Available methods

RungeKutta currently supports the following methods:

**Explicit**: `Euler`/`ExplicitEuler`, `Midpoint`/`ExplicitMidpoint`, `Heun2`, `Ralston2`, `Heun3`, `Kutta3`, `Ralston3`, `SSPRK3`, `RK4`, `Rule38`, `HeunEuler`, `Fehlberg45`/`F45`, `DormandPrince54`/`DP54`, `Verner65`/`V65`.

**Implicit**: `BackwardEuler`/`ImplicitEuler`, `ImplicitMidpoint`, `CrankNicolson`, `SDIRK3`, `GaussLegendre4`/`GL4`, `GaussLegendre6`/`GL6`, `LobattoIIIA4`, `LobattoIIIB2`, `LobattoIIIB4`, `LobattoIIIC2`, `LobattoIIIC4`, `RadauIA3`, `RadauIA5`, `RadauIIA3`, `RadauIIA5`.

<!-- <details><summary>Explicit</summary>

- `Euler`/`ExplicitEuler`
- `Midpoint`/`ExplicitMidpoint`
- `Heun2`
- `Ralston2`
- `Heun3`
- `Kutta3`
- `Ralston3`
- `SSPRK3`
- `RK4`
- `Rule38`
- `HeunEuler`
- `Fehlberg45`/`F45`
- `DormandPrince54`/`DP54`
- `Verner65`/`V65`

</details>

<details><summary>Implicit</summary>

- `BackwardEuler`/`ImplicitEuler`
- `ImplicitMidpoint`
- `CrankNicolson`
- `SDIRK3`
- `GaussLegendre4`/`GL4`
- `GaussLegendre6`/`GL6`
- `LobattoIIIA4`
- `LobattoIIIB2`
- `LobattoIIIB4`
- `LobattoIIIC2`
- `LobattoIIIC4`
- `RadauIA3`
- `RadauIA5`
- `RadauIIA3`
- `RadauIIA5`

</details> -->

## What's next?

Current plans for future developments are:

- Improve performance and error messages.
- Automatic size detection of stability region.
- IMEX methods, etc.
