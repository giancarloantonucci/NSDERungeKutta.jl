# NSDERungeKutta.jl

Serial Runge–Kutta integration: explicit, implicit, IMEX, and exponential methods.

## First example

```julia
using NSDERungeKutta
problem = IVP((u,t) -> [u[2], -sin(u[1])], [π/4, 0.0], (0.0, 2.0))
solver = RK4(h=0.01)
solution = solve(problem, solver)
@assert solution.t[end] == 2.0
solution(0.75)  # linear interpolation
```

Start with the [guided tour](guide.md), then consult:

- [Method catalogue and contracts](solvers.md)
- [Additional examples](examples.md)
- [API](api.md)

## Environment

These pages describe the current development checkout (Julia 1.10 or later).
Use compatible versions of the NSDE packages in one project. For local
checkouts, develop each package with `Pkg.develop(path="/path/to/package")`;
`Pkg.status()` reports what will actually load. A registered release may expose
an older API. Plotting examples additionally need Plots.
