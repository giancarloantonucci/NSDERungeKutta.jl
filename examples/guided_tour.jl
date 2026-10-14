# Executable companion to docs/src/guide.md.
# Run in an environment containing this package and the dependencies used below.

using NSDERungeKutta, LinearAlgebra
pendulum = IVP((u, t) -> [u[2], -sin(u[1])], [π/4, 0.0], (0.0, 2.0))
fixed = RK4(; h=0.03)
sol = solve(pendulum, fixed)
@assert sol.t[end] == 2.0
midpoint_value = sol(0.75)
@assert length(midpoint_value) == 2
nothing

adaptive = DormandPrince54(; h=0.02, εᵣ=1e-8, εₐ=1e-12, save_stepsizes=true)
a = solve(pendulum, adaptive)
@assert a.t[end] == 2.0
@assert length(a.stepsizes.accepted) == length(a.t)-1
@assert length(a.stepsizes.rejected) == length(a.stepsizes.accepted)+1
nothing

dense_sol = solve(pendulum, adaptive; dense=true)
continuous_value = dense_sol(0.751, adaptive.tableau)
@assert all(isfinite, continuous_value)
nothing

stiff = IVP([-1.0 0.0; 0.0 -100.0], [1.0, 1.0], (0.0, 1.0))
implicit = solve(stiff, RadauIIA5(; h=0.05))
@assert norm(implicit.u[end] - [exp(-1), exp(-100)]) < 1e-5

# This backward-Euler stage has no real solution: x = 1 + x^2.
impossible_stage = IVP((u, t) -> u .^ 2, [1.0], (0.0, 1.0))
failed = try
    solve(impossible_stage, BackwardEuler(; h=1.0, Mₙ=2))
    false
catch err
    err isa SolverFailure || rethrow()
    @assert err isa NewtonFailure
    true
end
@assert failed
nothing

using SparseArrays
n = 32
dx = 2π/n
x = dx .* collect(0:(n - 1))
D2 = spdiagm(-1=>ones(n-1), 0=>fill(-2.0, n), 1=>ones(n-1))
D2[1, n] = D2[n, 1] = 1.0
L = (0.01/dx^2) * D2
pde = IVP(L, (u, t) -> u .- u .^ 3, 0.1 .* sin.(x) .+ 0.05 .* cos.(3 .* x), (0.0, 0.2))
imex = solve(pde, IMEXSSP2_222(; h=0.01))
etd = solve(pde, ETDRK4(; h=0.01))
@assert all(isfinite, imex.u[end]) && all(isfinite, etd.u[end])
@assert norm(imex.u[end] - etd.u[end]) < 1e-3
nothing

cache = initialize_cache(pendulum, fixed)
reused = initialize_solution(pendulum, fixed)
for angle in (0.2, 0.4, 0.6)
    p = remake(pendulum; u0=[angle, 0.0])
    solve!(cache, reused, p, fixed)
    @assert reused.u[1] == p.u0
end
nothing
