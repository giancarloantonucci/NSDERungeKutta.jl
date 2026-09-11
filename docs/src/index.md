# NSDERungeKutta.jl

This is the documentation of [NSDERungeKutta.jl](https://github.com/giancarloantonucci/NSDERungeKutta.jl), a Julia package implementing Runge-Kutta methods: explicit, embedded adaptive explicit, diagonally implicit, fully implicit and implicit-explicit (IMEX) families over one shared solving shell.

## Installation

From the Julia REPL,

```
]add https://github.com/giancarloantonucci/NSDERungeKutta.jl
```

## Getting started

```julia
using NSDERungeKutta

problem = IVP((u, t) -> [u[2]; -sin(u[1])], [0.0, π/4], (0.0, 10.0))
solution = solve(problem, RK4(h = 1e-2))              # fixed step
solution = solve(problem, DormandPrince54(εᵣ = 1e-8)) # adaptive
u_mid = solution(5.0)                                  # interpolate
```

- The [Solvers](solvers.md) page lists the whole zoo with orders and explains adaptive stepping, dense output and the stability functions.
- The [Examples](examples.md) page walks through the bundled test problems with plots.
- The [API](api.md) holds the full reference.
