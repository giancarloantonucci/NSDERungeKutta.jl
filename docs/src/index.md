# RungeKutta.jl

This is the documentation of [RungeKutta.jl](https://github.com/giancarloantonucci/RungeKutta.jl), a Julia package implementing Runge-Kutta methods.

## Manual

```@contents
Depth = 3
```

## API

All exported types and functions are considered part of the public API and thus documented in this manual.

### Abstract types

```@autodocs
Modules = [RungeKutta]
Pages = ["abstract.jl"]
```

### Composite Types

```@autodocs
Modules = [RungeKutta]
Pages = [
  "tableau.jl",
  "stepsize.jl",
  "adaptive.jl",
  "newton.jl",
  "explicit.jl",
  "implicit.jl",
  "exponential.jl",
  "solution.jl",
  "cache.jl",
]
```

### Functions

```@docs
solve
solve!
step!
```

### Solvers

```@autodocs
Modules = [RungeKutta]
Pages = ["methods.jl"]
```

### Utilities

```@autodocs
Modules = [RungeKutta]
Pages = ["stability.jl"]
```

## Index

```@index
```
