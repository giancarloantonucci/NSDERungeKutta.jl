# RungeKutta

A Julia package implementing Runge-Kutta methods.

[![Build Status](https://github.com/antonuccig/RungeKutta.jl/workflows/CI/badge.svg)](https://github.com/antonuccig/RungeKutta.jl/actions)
[![Coverage](https://codecov.io/gh/antonuccig/RungeKutta.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/antonuccig/RungeKutta.jl)

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

We can plot the obtained `solution` by extracting its fields `u` and `t`, e.g. with the convenient macro `@↓ u, t = solution` from `ArrowMacros`. Alternatively, we can use the predefined recipes:

```julia
using Plots, LaTeXStrings
default(fontfamily = "Computer Modern")
plot(
  size = (800, 400),
  plot(solution, xlabel = L"t", label = [L"\theta" L"\omega"], legend = true),
  phaseplot(solution, vars = (1, 2), xlabel = L"\theta", ylabel = L"\omega")
)
# savefig("pendulum.svg")
```

![svg](images/pendulum.svg)

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

`RungeKutta` has also predefined recipes to plot stability regions and order stars:

```julia
plot(
  size = (1000, 400),
  stabilityf(RK4(), xlabel = L"\Re(z)", ylabel = L"\Im(z)", colour = :blues),
  orderstarf(RK4(), xlabel = L"\Re(z)", ylabel = L"\Im(z)", colour = :blues)
)
# savefig("regions.svg")
```

![svg](images/regions.svg)

<details><summary><b>Tests</b></summary>

[`OrdinaryDiffEq`](https://github.com/SciML/OrdinaryDiffEq.jl) is almost always faster and more efficient than `RungeKutta`, but `RungeKutta` is still competitive in a few cases, thanks to its simple design.

```julia
u0 = [2.0, 3.0, -14.0]
tspan = (0.0, 1.0)
problem = Lorenz(u0 = u0, tspan = tspan)
solver = RK4(h = 1e-3)
# solver = DP54(h = 1e-3)
# solver = CrankNicolson(h = 1e-3)
```

```julia
using BenchmarkTools
@benchmark solve($problem, $solver)
```

```julia
BenchmarkTools.Trial:
  memory estimate:  313.64 KiB
  allocs estimate:  6000
  --------------
  minimum time:     315.252 μs (0.00% GC)
  median time:      337.454 μs (0.00% GC)
  mean time:        470.207 μs (7.79% GC)
  maximum time:     17.859 ms (0.00% GC)
  --------------
  samples:          10000
  evals/sample:     1
```

```julia
# import OrdinaryDiffEq; const DE = OrdinaryDiffEq # v1.0
import OrdinaryDiffEq as DE #v1.6
u0 = [2.0, 3.0, -14.0]
tspan = (0.0, 1.0)
function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - 8/3 * u[3]
end
DE_problem = DE.ODEProblem(lorenz!, u0, tspan)
DE_solver = DE.RK4()
# DE_solver = DE.DP5()
# DE_solver = DE.Trapezoid()
```

```julia
@benchmark DE.solve($DE_problem, $DE_solver, dt = 1e-3, adaptive = false)
```

```julia
BenchmarkTools.Trial:
  memory estimate:  448.80 KiB
  allocs estimate:  4034
  --------------
  minimum time:     369.874 μs (0.00% GC)
  median time:      394.468 μs (0.00% GC)
  mean time:        489.820 μs (9.59% GC)
  maximum time:     12.049 ms (95.71% GC)
  --------------
  samples:          9973
  evals/sample:     1
```

</details>

## Methods

`RungeKutta` currently supports the following methods:

<!-- explicit (`Euler`/`ExplicitEuler`, `Midpoint`/`ExplicitMidpoint`, `Heun2`, `Ralston2`, `Heun3`, `Kutta3`, `Ralston3`, `SSPRK3`, `RK4`, `Rule38`, `HeunEuler`, `Fehlberg45`/`F45`, `DormandPrince54`/`DP54`, `Verner65`/`V65`) and implicit methods (`BackwardEuler`/`ImplicitEuler`, `ImplicitMidpoint`, `CrankNicolson`, `SDIRK3`, `GaussLegendre4`/`GL4`, `GaussLegendre6`/`GL6`, `LobattoIIIA4`, `LobattoIIIB2`, `LobattoIIIB4`, `LobattoIIIC2`, `LobattoIIIC4`, `RadauIA3`, `RadauIA5`, `RadauIIA3`, `RadauIIA5`). -->

<details><summary>Explicit</summary>

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

</details>

## What's next?

Current plans for future developments are:
- Improve performance and error messages.
- Automatic size detection of stability region.
- IMEX methods.
