# Changelog

## 0.2.0

Requires NSDEBase 0.3.1 (four-argument RHS form).

### Fixed (numerical)
- Butcher tableaus: `Butcher7` coefficient `5/164` → `5/154`; `RadauIIA5`
  denominators `75` → `225` in two entries; `IMEXSSP3_332` implicit tableau
  rows and weights corrected, and its order stated as 2 (SSP3 names the
  explicit part).
- Adaptive stepping: the retry counter was reset before the budget check, so
  the "maximum reductions" warning could never fire; `StepSizes.accepted`
  recorded the controller's next proposal instead of the step taken
  (`hs.accepted == diff(solution.t)` now holds); a zero error estimate no
  longer divides by zero.
- Dense output: stage history is stored per accepted step and sized to the
  intervals, so `sol[end]` on a dense solution no longer throws.
- A fixed-step schedule that divides the span ends on `tN` bit for bit
  (`fixedstep_count`); schedules that do not divide it keep the overshoot
  policy; adaptive solvers are untouched.
- `RK3` documented with its own docstring (it showed `RungeKutta4`'s).
- Plots recipe skips the trailing empty rejection bin.

### Changed (breaking)
- Implicit solvers (DIRK, IRK, IERK) accept a Newton stage on its residual,
  `‖r‖ ≤ εₐ + εᵣ·max(‖x‖, ‖f‖)`, with `hypot`-accumulated norms; non-finite
  iterates are rejected; exhausting `Mₙ` throws `NewtonFailure`. No
  unconverged stage is ever used. `NewtonParameters` gains `εₐ` (default
  `1e-12`), validates its inputs, and its default `εᵣ` moves from `1e-3` to
  `1e-8` — an accuracy choice for a bound on the iterate actually used; cost
  not benchmarked. All named implicit solvers accept `εₐ`. IERK's stiff stage
  value is evaluated at the accepted iterate.
- Dense `RungeKuttaSolution` slicing follows the interval storage: a node has
  no stages, a range has `length − 1`, non-contiguous slices throw; writes
  into a dense solution through the indexer must replace the whole solution.
- One family-agnostic `adaptivestep!`; `AdaptiveParameters` on a tableau
  without an embedded pair throws at construction; IMEX solvers refuse
  adaptive parameters.
- `solve(problem, solver; dense)` no longer forwards arbitrary keywords.
- Exponential Runge–Kutta rewritten: `ExponentialRungeKuttaSolver`/`EXPRK`
  on `ExponentialTableau` (operator coefficients, EXPINT convention) replaces
  `ExplicitExponentialRungeKuttaSolver`/`ExRK`. `ETDRK4` keeps its name.
- Supported Julia: `1.6` and later (was `1.10`). Verified on 1.6–1.13. On
  1.6–1.8 the CHOLMOD factor method of `directldiv!` is taken from the
  `SuiteSparse` stdlib in the system image; from 1.9 it is
  `SparseArrays.CHOLMOD`.

### Added
- Exponential solvers `LawsonEuler`, `Lawson4`, `NorsettEuler`/`ETDEuler`/
  `ExponentialEuler`, `ETD2RK`, `ETD3RK`, `ETD4RK`, `Krogstad`,
  `HochbruckOstermann4`/`HochOst4`; `phifunctions`, `expphifunctions`.
- `NewtonParameters` and `NewtonFailure` exported.
- Tableau pretty-printing; `directldiv!` for CHOLMOD factors, which lack an
  in-place `ldiv!`.
- Docs: solver guide, API page; Aqua in the test suite; the exponential-RK
  oracle runs under `Pkg.test`.

### Removed
- The unused `MakieCore` and `Pkg` dependencies.

### Known limits
- A Newton failure inside an adaptive implicit step is thrown before the
  controller can reject the step; automatic retry is not implemented.

### Migration
- Replace `ExRK`/`ExplicitExponentialRungeKuttaSolver` with `EXPRK`/
  `ExponentialRungeKuttaSolver`.
- Code that relied on an implicit solver silently returning after `Mₙ`
  iterations will now see `NewtonFailure`; shorten `h`, raise `Mₙ`, or loosen
  the tolerances deliberately.
- `sol[i]` on a dense solution no longer carries a stage entry; slice a range.
