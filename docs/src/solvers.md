# Solvers

One shared shell — tableau, step size, the solve loop, dense output, stability functions — with each family owning only its `step!` kernel, dispatched on the right-hand-side type.

## Explicit (ERK)

`Euler`/`ExplicitEuler` (1), `Heun2`, `Midpoint`/`ExplicitMidpoint`, `Ralston2` (2), `Heun3`, `RungeKutta3`/`RK3`, `Ralston3`, `SSPRK3` (3), `Ralston4`, `RungeKutta4`/`RK4`, `Rule38` (4), `Butcher5`, `KuttaNystrom5` (5), `Butcher6` (6), `Butcher7` (7).

## Embedded explicit, adaptive

`HeunEuler` 2(1), `BogackiShampine` 3(2), `Fehlberg45`, `DormandPrince54`, `Verner65`, `Fehlberg78`.

Adaptive stepping needs an embedded pair, so it is available for these solvers only: pass tolerances (`εᵣ`, `εₐ`) and, if wanted, `save_stepsizes = true` to record accepted and rejected steps. Handing `AdaptiveParameters` to a solver whose tableau has no embedding throws at construction — no silent fixed stepping.

## Diagonally implicit (DIRK)

`BackwardEuler`/`ImplicitEuler` (1), `ImplicitMidpoint`/`GaussLegendre2`, `SDIRK2`, `LobattoIII2`, `CrankNicolson`/`LobattoIIIA2` (2), `SDIRK3`, `RadauI3`, `RadauII3` (3), `SDIRK4`, `LobattoIII4` (4).

Simplified Newton per stage; a direct linear solve when the right-hand side is a `LinearRightHandSide`.

Every implicit family shares one Newton contract, set through `εᵣ`, `εₐ` and `Mₙ` (see `NewtonParameters`): a stage is accepted when its **residual** — how far the current iterate is from solving the stage equation — satisfies `‖r‖ ≤ εₐ + εᵣ·max(‖x‖, ‖f‖)`. The size of the last update is never the verdict, because a tiny update is also what a stalled or badly scaled iteration produces. Non-finite iterates and residuals are rejected outright, tolerances must be finite, and if `Mₙ` updates pass without acceptance the step **throws a `NewtonFailure`** naming the time, the stage, the residual and the tolerance it missed. An unconverged stage is never used silently: a step whose stage equation has no solution (say backward Euler with `h = 1` on `u′ = u²`, `u₀ = 1`) fails loudly rather than returning a finite number with a large residual. Catch the exception to shorten the step or raise `Mₙ`. One limit to know about: when a generic `DIRK` is built with an embedded tableau and `AdaptiveParameters`, a Newton failure is thrown before the controller can reject and shorten the step — automatic recovery for adaptive implicit stepping is not part of this release. The named implicit solvers are fixed-step and unaffected.

## Fully implicit (IRK)

`LobattoIIIC2` (2), `RadauIA3`, `RadauIIA3` (3), `GaussLegendre4`, `LobattoIIIA4`, `LobattoIIIB4`, `LobattoIIIC4` (4), `RadauI5`, `RadauIA5`, `RadauII5`, `RadauIIA5` (5), `GaussLegendre6` (6).

One coupled Newton solve over all stages per step.

## Implicit-explicit (IMEX)

`IMEXEuler` (1), `IMEXSSP2_222`, `IMEXSSP2_322`, `IMEXSSP2_332` (2), `IMEXSSP3_332` (2 — the SSP**3** in the name refers to the SSP order of the explicit part).

These solve a `SplitRightHandSide` ``f = f_s + f_{ns}``: the stiff part implicitly (direct solve if linear, Newton if not), the non-stiff part explicitly.

## Dense output and interpolation

`solution(t)` interpolates: linear spline by default, cubic Hermite given the derivative (`solution(t, f)`), and — with `solve(…; dense = true)` on a solver whose tableau carries dense-output weights — the method's own continuous extension via `solution(t, solver.tableau)`.

## Stability functions

`stability_function(z, solver)` evaluates ``R(z)`` (scalar or matrix argument); the plot recipes `stability`/`stabilityf` and `orderstar`/`orderstarf` draw stability regions and order stars.
