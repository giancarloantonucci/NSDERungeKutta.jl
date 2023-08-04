# NSDERungeKutta.jl

A Julia package implementing Runge-Kutta methods.

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://giancarloantonucci.github.io/NSDERungeKutta.jl/dev) ![Build Status](https://img.shields.io/github/actions/workflow/status/giancarloantonucci/NSDERungeKutta.jl/CI.yml) ![Coverage Status](https://img.shields.io/codecov/c/github/giancarloantonucci/NSDERungeKutta.jl)

## Installation

This package is a [registered package](https://juliahub.com/ui/Search?q=NSDERungeKutta&type=packages) compatible with Julia v1.6 and above. From the Julia REPL,

```
]add NSDERungeKutta
```

Read the [documentation](https://giancarloantonucci.github.io/NSDERungeKutta.jl/dev) for a complete overview of this package.

## Usage

Let's say that we want to solve the [simple gravity pendulum problem](https://en.wikipedia.org/wiki/Pendulum_(mathematics)#Simple_gravity_pendulum) using the [midpoint method](https://en.wikipedia.org/wiki/Midpoint_method). Here is one way to do it with NSDERungeKutta.jl:

```julia
using NSDERungeKutta
f(u, t) = [u[2]; -sin(u[1])]
u0 = [0.0; π/4]
tspan = (0.0, 2π/3 * √9.81)
problem = IVP(f, u0, tspan)
solver = Midpoint(h = 1e-2)
solution = solve(problem, solver)
```

We can plot the obtained `solution` by extracting its fields `u` and `t`, e.g. using the convenient macro `@↓ u, t = solution` from [ArrowMacros.jl](https://github.com/giancarloantonucci/ArrowMacros.jl). Alternatively, we can use the available predefined recipes:

```julia
using Plots, LaTeXStrings
gr(fontfamily = "Computer Modern", framestyle = :box, label = "", tickdirection = :out)
p₁ = plot(solution, xlabel = L"t", label = [L"\theta" L"\omega"])
p₂ = phaseplot(solution, variables = (1, 2), xlabel = L"\theta", ylabel = L"\omega")
plot(size = (900, 450), p₁, p₂, left_margin = 3Plots.mm, bottom_margin = 3Plots.mm)
```

![svg](imgs/pendulum.svg)

For convenience, this package re-exports all the ODE problems defined in [NSDEBase.jl](https://github.com/giancarloantonucci/NSDEBase.jl), e.g. `SimplePendulum` for the above problem.

This package has some predefined recipes to plot **stability regions** and **order stars** too:

```julia
p₁ = stabilityf(RK4(), xlabel = L"\Re(z)", ylabel = L"\Im(z)", colour = :blues, resolution = 500)
p₂ = orderstarf(RK4(), xlabel = L"\Re(z)", ylabel = L"\Im(z)", colour = :blues, resolution = 500)
plot(size = (1000, 400), p₁, p₂, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
```

![svg](imgs/regions.svg)

## Methods

This package currently supports the following methods:

**Explicit**:

- `Euler`/`ExplicitEuler`, `Heun2`, `Midpoint`/`ExplicitMidpoint`, `Ralston2`, `Heun3`, `RungeKutta3`/`RK3`, `Ralston3`, `SSPRK3`, `Ralston4`, `RungeKutta4`/`RK4`, `Rule38`, `Butcher5`, `KuttaNystrom5`, `Butcher6`, `Butcher7`,
- (Embedded) `HeunEuler`, `BogackiShampine`, `Fehlberg45`, `DormandPrince54`, `Verner65`, `Fehlberg78`.

**Diagonally Implicit**:
- `BackwardEuler`/`ImplicitEuler`, `ImplicitMidpoint`/`GaussLegendre2`, `SDIRK2`, `LobattoIII2`, `CrankNicolson`/`LobattoIIIA2`, `SDIRK3`, `RadauI3`, `RadauII3`, `SDIRK4`, `LobattoIII4`.

**Implicit**:
- `LobattoIIIC2`, `RadauIA3`, `RadauIIA3`, `GaussLegendre4`, `LobattoIIIA4`, `LobattoIIIB4`, `LobattoIIIC4`, `RadauI5`, `RadauIA5`, `RadauII5`, `RadauIIA5`, `GaussLegendre6`.
