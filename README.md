# NSDERungeKutta.jl

Serial Runge–Kutta integration: explicit, implicit, IMEX, and exponential methods.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://giancarloantonucci.github.io/NSDERungeKutta.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://giancarloantonucci.github.io/NSDERungeKutta.jl/dev) ![Build Status](https://img.shields.io/github/actions/workflow/status/giancarloantonucci/NSDERungeKutta.jl/CI.yml) ![Coverage Status](https://img.shields.io/codecov/c/github/giancarloantonucci/NSDERungeKutta.jl)

## Start here

```julia
using NSDERungeKutta
problem = IVP((u, t) -> [u[2], -sin(u[1])], [π/4, 0.0], (0.0, 2.0))
solver = RK4(h=0.01)
solution = solve(problem, solver)
@assert solution.t[end] == 2.0
solution(0.75)  # linear interpolation
```

All named integration methods require `h`. For an adaptive method such as
`DormandPrince54(h=0.01, εᵣ=1e-8, εₐ=1e-12)`, it is the first trial step.
Fixed solves shorten the last step when needed. `solution(t)` uses linear
interpolation; optional method-specific continuous output requires `dense=true`
and `solution(t, solver.tableau)`. See the guide before enabling it.

NSDEBase is re-exported. `NewtonFailure` and `StepFailure` identify numerical
failures. Reuse `initialize_cache`, `initialize_solution`, and
`solve!(cache, solution, problem, solver)` to reduce setup work.

See [Solvers](docs/src/solvers.md) for the method catalogue, endpoint rules,
tolerances, dense-output support, and exponential-operator costs. The
[additional examples](docs/src/examples.md) cover built-in problems.

## Learn by running

The [guided tour](docs/src/guide.md) explains each step. Its companion
[`examples/guided_tour.jl`](examples/guided_tour.jl) runs the numerical examples
with assertions. Run it in a project containing this package and its example
dependencies; use Julia 1.10 or later.

The ecosystem is NSDEBase → NSDERungeKutta → NSDETimeParallel → NSDEMovingWindow.
ArrowMacros provides supporting macros; they are not required in user examples.

## Installation and versions

For a registered release, use `Pkg.add("NSDERungeKutta")`. These docs describe the
current development API; a registered version may differ. To use a local
checkout, activate your working environment and `Pkg.develop(path="/path/to/NSDERungeKutta")`.
Develop matching local dependency checkouts in that environment too, then check
`Pkg.status()`. Plotting requires Plots.

See the [documentation](https://giancarloantonucci.github.io/NSDERungeKutta.jl/dev/)
and [API reference](docs/src/api.md) for details.

## Build the manual locally

From the package root, with matching dependency checkouts already available:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs --threads=4 docs/make.jl
```

Documenter executes the guided-tour examples during the build. For local
unreleased dependencies, develop their paths into `docs` before instantiating.
The deployment step only publishes in an appropriately configured CI environment.
