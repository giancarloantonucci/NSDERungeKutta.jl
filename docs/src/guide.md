# Guided tour: serial integration

This page moves from a small non-stiff problem to a split spatial discretisation.
The executable companion is `examples/guided_tour.jl`. These examples check API
usage and basic numerical behaviour; they are not performance benchmarks.

## Fixed steps and interpolation

```@example serial
using NSDERungeKutta, LinearAlgebra
pendulum = IVP((u,t) -> [u[2], -sin(u[1])], [π/4, 0.0], (0.0, 2.0))
fixed = RK4(h=0.03)
sol = solve(pendulum, fixed)
@assert sol.t[end] == 2.0
midpoint_value = sol(0.75)
@assert length(midpoint_value) == 2
nothing # hide
```

`sol.t` and `sol.u` hold the stored nodes and states. `sol(t)` uses linear
interpolation, regardless of the order of the time integrator. The last step is
shortened to reach the endpoint when necessary; negligible remainders are folded
into the preceding step. See [Solvers](solvers.md) for the precise endpoint rule.

## Adaptive steps

```@example serial
adaptive = DormandPrince54(h=0.02, εᵣ=1e-8, εₐ=1e-12, save_stepsizes=true)
a = solve(pendulum, adaptive)
@assert a.t[end] == 2.0
@assert length(a.stepsizes.accepted) == length(a.t)-1
@assert length(a.stepsizes.rejected) == length(a.stepsizes.accepted)+1
nothing # hide
```

`h` is the first trial step and is required. The embedded local error is
`h*sum((b-d)*k)`. Each component is divided by
`εₐ + εᵣ*max(abs(u_old),abs(u_new))`; the RMS scaled error must be at most one.
This controls an estimate of local error, not the global trajectory error.
`rejected[j]` contains the rejected trial sizes before accepted step `j`.
The history also keeps one trailing empty retry bin for the next step.

At the retry limit, the current implementation can accept a **finite** step
outside tolerance with a warning. A nonfinite trial is shrunk and retried;
exhausted nonfinite retries or inability to advance time raise `StepFailure`.
Do not discard those warnings when judging accuracy.

## Optional continuous output

```@example serial
dense_sol = solve(pendulum, adaptive; dense=true)
continuous_value = dense_sol(0.751, adaptive.tableau)
@assert all(isfinite, continuous_value)
nothing # hide
```

The method's continuous extension needs stored stages and the tableau argument.
A bare `dense_sol(t)` still uses linear interpolation. Only Euler (order 1),
BogackiShampine (order 3), and DormandPrince54/DP54 (order 4) provide these
extensions. Other solvers reject `dense=true`. It is optional and uses more
memory; ordinary integration does not need it.

## Stiffness and numerical failures

```@example serial
stiff = IVP([-1.0 0.0; 0.0 -100.0], [1.0, 1.0], (0.0, 1.0))
implicit = solve(stiff, RadauIIA5(h=0.05))
@assert norm(implicit.u[end] - [exp(-1), exp(-100)]) < 1e-5

# This backward-Euler stage has no real solution: x = 1 + x^2.
impossible_stage = IVP((u,t) -> u.^2, [1.0], (0.0, 1.0))
failed = try
    solve(impossible_stage, BackwardEuler(h=1.0, Mₙ=2))
    false
catch err
    err isa SolverFailure || rethrow()
    @assert err isa NewtonFailure
    true
end
@assert failed
nothing # hide
```

An implicit stability region does not guarantee that nonlinear stage equations
can be solved at an arbitrary step. Inspect `NewtonFailure`, reduce `h`, or
change the method; increasing `Mₙ` alone need not help. `NewtonFailure` and
`StepFailure` are both numerical `SolverFailure`s. Catch only the failures you
intend to handle and rethrow unrelated exceptions.

## A split PDE: periodic Allen–Cahn

Discretise ``u_t=εu_{xx}+u-u^3`` on a periodic grid. The matrix below represents
diffusion, and the function represents the reaction.

```@example serial
using SparseArrays
n = 32
dx = 2π/n
x = dx .* collect(0:n-1)
D2 = spdiagm(-1=>ones(n-1), 0=>fill(-2.0,n), 1=>ones(n-1))
D2[1,n] = D2[n,1] = 1.0
L = (0.01/dx^2) * D2
pde = IVP(L, (u,t) -> u .- u.^3,
          0.1 .* sin.(x) .+ 0.05 .* cos.(3 .* x), (0.0, 0.2))
imex = solve(pde, IMEXSSP2_222(h=0.01))
etd = solve(pde, ETDRK4(h=0.01))
@assert all(isfinite, imex.u[end]) && all(isfinite, etd.u[end])
@assert norm(imex.u[end] - etd.u[end]) < 1e-3
nothing # hide
```

The diffusion eigenvalues lie in `[-4ε/dx^2,0]`; refining the spatial grid
increases stiffness. IMEX handles diffusion implicitly, while exponential
methods use functions of the linear operator. Neither removes all restrictions
from the explicitly treated nonlinear term. Refine both space and time before
interpreting a PDE result.

The current exponential implementation has an elementwise path for diagonal
operators. General sparse operators are converted to dense matrices for operator
functions. This small example is suitable for learning; it is not a matrix-free
large-PDE solver or a scalability claim.

## Reuse storage and measure the right work

```@example serial
cache = initialize_cache(pendulum, fixed)
reused = initialize_solution(pendulum, fixed)
for angle in (0.2, 0.4, 0.6)
    p = remake(pendulum; u0=[angle, 0.0])
    solve!(cache, reused, p, fixed)
    @assert reused.u[1] == p.u0
end
nothing # hide
```

The result buffer is overwritten on each solve. Copy any final state you want to
retain. Warm up before timing, and report whether initialization is included.
Cache reuse reduces setup work but does not ensure an allocation-free solve.
Rebuild operator-dependent caches when the linear operator or method changes.

## Plotting

With Plots installed, these commands use the variables defined above:

```julia
using Plots
plot(sol; variables=[1,2], label=["angle" "velocity"])
phaseplot(sol; variables=(1,2))
plot(a.stepsizes)
stabilityf(RK4)
orderstarf(RadauIIA5)
```

Stability recipes accept ordinary RK constructors, configured solvers, tableaus,
or scalar stability functions. IMEX and exponential solver objects need a
problem-specific stability analysis and are refused by these scalar recipes.
