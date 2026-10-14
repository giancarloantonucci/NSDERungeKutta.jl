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

## Exponential (EXPRK)

`LawsonEuler`, `NorsettEuler`/`ETDEuler`/`ExponentialEuler` (1), `ETD2RK` (2), `ETD3RK` (3), `ETD4RK`/`ETDRK4`, `Lawson4`, `Krogstad`, `HochbruckOstermann4`/`HochOst4` (4).

For semilinear problems ``u' = Lu + g(t) + f_{ns}(u, t)``, supplied as a `SplitRightHandSide` whose stiff part is a `LinearRightHandSide` (or as a plain `LinearRightHandSide`; exact propagation of the homogeneous linear part does not imply exact integration of arbitrary time-dependent forcing). The linear part is propagated exactly: the stage and weight coefficients of an `ExponentialTableau` are operator functions of ``z = hL`` built from the φ-functions ``\varphi_k``, which `phifunctions` evaluates by scaled diagonal Padé (the EXPINT convention). Lawson schemes are classical Runge–Kutta methods in the integrating-factor variables; ETD, Krogstad and Hochbruck–Ostermann schemes are exponential integrators proper, with the stiff-order conditions of the latter. Because ``e^{hL}`` and the ``\varphi_k(hL)`` are computed in each solver cache, the step size is fixed: these solvers take no adaptive parameters.

The tests in `test/exponential_oracles.jl` and `test/exponential_methods.jl` check the φ engine against 256-bit ground truth and a contour integral, the shipped tableaus against their operator identities, and the whole against the Kuramoto–Sivashinsky benchmark.

## Dense output and interpolation

IMEX solvers require `dense=false`: storing and interpolating both implicit
and explicit stage families is not yet implemented. A request for `dense=true`
throws `ArgumentError` before stepping; ordinary linear interpolation remains
available on the returned solution.

`solution(t)` interpolates: linear spline by default, cubic Hermite given the derivative (`solution(t, f)`), and — with `solve(…; dense = true)` on a solver whose tableau carries dense-output weights — the method's own continuous extension via `solution(t, solver.tableau)`.

## Stability functions

`stability_function(z, solver)` evaluates ``R(z)`` (scalar or matrix argument); the plot recipes `stability`/`stabilityf` and `orderstar`/`orderstarf` draw stability regions and order stars.

## Endpoint and run-state contract

Every integration method requires an explicit positive step size, e.g.
`RK4(h=0.01)` or `DormandPrince54(h=0.01, εᵣ=1e-8)`. Stability and order-star
plots accept a method constructor (`stabilityf(RK4)`), a tableau, a configured
solver, or a scalar stability function; no integration step is needed for them.

Fixed-step solves shorten the final step to land on the endpoint. A remainder
no larger than `1e-6*h` is folded into the preceding step. Adaptive solves use
the same small-remainder rule, including a floating-point endpoint allowance.
A step that cannot advance time raises a `StepFailure`: the requested time
resolution is then unavailable in that floating-point type.

`chunk_solver(solver, duration)` prepares a nominal fixed schedule. Parareal
calls this automatically; window shifts retain the nominal duration until the
window length or solver changes. Each chunk cache retains any exponential
coefficients for both the full and final step. Coefficients must be rebuilt
when the linear operator, step size, or nominal duration changes.

Solver settings are immutable. `withstepsize(solver, h)` returns changed
settings; adaptive integration passes each realised step directly to the
kernel. Step histories belong to `solution.stepsizes`, and all mutable stage
buffers belong to a particular solve.

## Numeric types and adaptive error control

Problem construction converts integer/rational states and time endpoints to floating
point and preserves already-floating state storage. The time grid and compensated
summation accumulator use the problem's time type, including Float32 and BigFloat;
the supplied step is converted to that type when advancing time. This preserves
storage precision; it does not upgrade Float64 tableau coefficients to arbitrary
precision. Choose a complex initial state for complex forcing, even with a real
linear operator.

Embedded methods estimate the local step error `h*sum((b-d)*k)`. Each component is
scaled by `εₐ + εᵣ*max(abs(u_old),abs(u_new))`; the RMS of these scaled errors is
accepted when it is at most one. Repeating identical components therefore leaves
the effective tolerance unchanged. A zero scale accepts only zero estimated error.
The default absolute tolerance is 1e-12; tolerances must be finite, nonnegative,
and not both zero. Step multipliers use safety 0.9 and are clamped to [0.2,2].
At the retry limit a finite rejected step is still accepted with a warning; it has
not met the requested tolerance. Nonfinite estimates cannot be accepted this way; exhausted nonfinite retries raise `StepFailure`.
Changing adaptive tolerances or controller settings changes step histories and numerical results.

Plain and scheduled exponential solves retain coefficients for the full step and
one final-step length. Reusing a cache over the same span reuses both. A changed
final-step length replaces the final-step cache. The cache remains tied to its
linear operator and main step; rebuild it when either changes.

### Continuous extensions and numerical checks

`dense=true` is supported by **Euler** (order 1), **BogackiShampine** (order 3)
and **DormandPrince54 / DP54** (order 4). Other methods reject this request
before solving; `dense=false` retains ordinary linear interpolation without
storing unused stages. Call `solution(t, solver.tableau)` to use the extension.
The Dormand–Prince coefficients are expanded from Hairer and Wanner's
[DOPRI5 continuous output](https://www.unige.ch/~hairer/prog/nonstiff/dopri5.f),
described in Hairer, Nørsett and Wanner, *Solving Ordinary Differential
Equations I*, second edition (1993). Tests check rooted-tree conditions at
interior fractions as well as the endpoint weights.

The test-only [RootedTrees.jl](https://sciml.github.io/RootedTrees.jl/stable/tutorials/RK_order_conditions/)
audit checks every primary and embedded ordinary RK tableau through its
claimed order and checks the next order fails. IMEX coupling and exponential
methods have separate tests; ordinary rooted trees do not certify those families.

Nonfinite adaptive trial estimates are rejected with the minimum step factor
0.2. Retry exhaustion or a step too small to advance time still throws.
Rejected trials restore the time-compensation accumulator. Implicit Newton
iteration starts from the stage predictor in the nonlinear IMEX kernel.
Newton tolerances cannot both be zero.

Scalar stability plots accept ordinary Butcher tableaus, corresponding method
constructors/solvers, or a user-supplied scalar stability function. They reject
IMEX and exponential solver objects with `ArgumentError`.
