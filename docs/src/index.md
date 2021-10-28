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

### Composite types

```@autodocs
Modules = [RungeKutta]
Pages = [
  "explicit/struct.jl",
  "implicit/struct.jl",
  "exponential/struct.jl",
  "solution.jl",
  "stepsize.jl",
  "tableau.jl",
  "newton.jl",
  "adaptive.jl",
  "explicit/cache.jl",
  "implicit/cache.jl",
  "exponential/cache.jl"]
```

#### Solvers

```@autodocs
Modules = [RungeKutta]
Pages = [
  "explicit/solvers.jl",
  "implicit/solvers.jl",
  "exponential/solvers.jl"]
```

### Functions

```@autodocs
Modules = [RungeKutta]
Pages = [
  "solve.jl",
  "explicit/step.jl",
  "implicit/step.jl",
  "exponential/step.jl"]
```

### Utilities

```@autodocs
Modules = [RungeKutta]
Pages = ["stability.jl"]
```

## Index

```@index
```
