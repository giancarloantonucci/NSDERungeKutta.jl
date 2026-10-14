# API: types and functions

The exported types and functions are the public API. The internal helpers are
listed for reference and may change without notice. The named methods are on
their own pages: [explicit](api_explicit.md), [implicit](api_implicit.md), and
[IMEX and exponential](api_imex_exponential.md).

```@autodocs
Modules = [NSDERungeKutta]
Pages = [
    "/src/NSDERungeKutta.jl",
    "/src/abstract.jl",
    "/src/adaptive.jl",
    "/dirk/cache.jl",
    "/dirk/constructor.jl",
    "/dirk/step.jl",
    "/erk/cache.jl",
    "/erk/constructor.jl",
    "/erk/step.jl",
    "/exprk/cache.jl",
    "/exprk/constructor.jl",
    "/exprk/step.jl",
    "/exprk/tableau.jl",
    "/ierk/cache.jl",
    "/ierk/constructor.jl",
    "/ierk/step.jl",
    "/irk/cache.jl",
    "/irk/constructor.jl",
    "/irk/step.jl",
    "/src/newton.jl",
    "/src/phi.jl",
    "/src/plots_recipes.jl",
    "/src/schedule.jl",
    "/src/solution.jl",
    "/src/solve.jl",
    "/src/stability.jl",
    "/src/stepsize.jl",
    "/src/tableau.jl",
    "/src/tableau_show.jl",
    "/src/utils.jl",
]
```
