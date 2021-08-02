# RungeKutta.jl

```@contents
```

## Public API

### Constructors

```@docs
solve
solve!
RungeKuttaSolver
ExplicitRungeKuttaSolver
ImplicitRungeKuttaSolver
RungeKuttaSolution
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
ButcherTableau
ℛ
```

## Index

```@index
```
