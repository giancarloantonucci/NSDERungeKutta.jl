# RungeKutta.jl

This is the documentation of [RungeKutta.jl](https://github.com/antonuccig/RungeKutta.jl), a Julia package implementing Runge-Kutta methods.

## Manual

```@contents
Depth = 3
```

## API

All exported types and functions are considered part of the public API and thus documented in this manual.

### Types

```@docs
ButcherTableau
RungeKuttaSolver
ExplicitRungeKuttaSolver
ImplicitRungeKuttaSolver
RungeKuttaSolution
```

### Functions

```@docs
solve
solve!
```

### Solvers

#### Explicit

```@docs
Euler
Midpoint
Heun2
Ralston2
Heun3
Kutta3
Ralston3
SSPRK3
RK4
Rule38
HeunEuler
Fehlberg45
DormandPrince54
Verner65
```

#### Implicit

```@docs
BackwardEuler
ImplicitMidpoint
CrankNicolson
SDIRK3
GaussLegendre4
GaussLegendre6
LobattoIIIA4
LobattoIIIB2
LobattoIIIB4
LobattoIIIC2
LobattoIIIC4
RadauIA3
RadauIA5
RadauIIA3
RadauIIA5
```

### Utilities

```@docs
ℛ
```

## Index

```@index
```
