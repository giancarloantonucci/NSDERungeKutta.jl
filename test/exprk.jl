# NSDERungeKutta/test/exprk.jl
#
# Include from the END of runtests.jl (`include("exprk.jl")`): this file uses
# the `convergence_order`/`sweep_params` harness and the `AP` const defined
# there. SparseArrays resolves through the package's own dependencies.

using SparseArrays: spdiagm

@testset "φ-functions (phipade port)" begin
    # Closed forms at scalar arguments: real, negative-large, complex.
    for z in (0.3, -1.7, 2.5 + 1.3im, -40.0)
        φ = phifunctions(z, 3)
        @test φ[1] ≈ expm1(z) / z                              rtol = 1e-13
        @test φ[2] ≈ (exp(z) - 1 - z) / z^2                    rtol = 1e-12
        @test φ[3] ≈ (exp(z) - 1 - z - z^2 / 2) / z^3          rtol = 1e-12
    end

    # Removable singularity: φₖ(0) = 1/k! comes straight out of the Padé
    # normalisation, no limit-taking involved.
    @test phifunctions(0.0, 4) ≈ [1.0, 1/2, 1/6, 1/24]

    # Matrix arguments: the recurrence z φₖ₊₁ = φₖ − I/k! and consistency with
    # the matrix exponential, on a deterministic dense matrix.
    M = [sin(3i + j) for i = 1:5, j = 1:5] - 3.0I
    φ = phifunctions(M, 4)
    for k = 1:3
        @test M * φ[k+1] ≈ φ[k] - Matrix(1.0I, 5, 5) / factorial(k)  atol = 1e-12
    end
    @test M * φ[1] + I ≈ exp(M)  rtol = 1e-12

    # The elementwise Diagonal fast path must agree with the dense path,
    # including through the scaling-and-squaring branch (|λ| = 100 forces it).
    λ = [-100.0, -1.0, 0.0, 2.0]
    φd = phifunctions(Diagonal(λ), 3)
    φm = phifunctions(Matrix(Diagonal(λ)), 3)
    for k = 1:3
        @test Matrix(φd[k]) ≈ φm[k]  atol = 1e-13
    end
end

@testset "chain exponential (expphifunctions) under strong damping" begin
    # The identity zφ₁(z) + I is only ABSOLUTELY accurate: at z = -40 it
    # returns 0 with unbounded relative error. The chain exponential must
    # stay relatively accurate down to the underflow threshold.
    for z in (-40.0, -200.0, -700.0)
        ez, = expphifunctions(z, 1)
        @test ez ≈ exp(z)  rtol = 1e-11
    end
    λ = [-240.0, -40.0, -1.0, 0.0]
    ezd, φd = expphifunctions(Diagonal(λ), 2)
    @test all(isapprox.(ezd.diag, exp.(λ); rtol = 1e-11)) # elementwise, incl. e⁻²⁴⁰ ≈ 3e-105
    # phifunctions is exactly the φ part of expphifunctions:
    @test phifunctions(-3.3, 3) == expphifunctions(-3.3, 3)[2]
    @test Matrix(φd[1]) ≈ phifunctions(Matrix(Diagonal(λ)), 1)[1]  atol = 1e-13
end

@testset "stiff modes decay to their true exponential (regression)" begin
    # With the old zφ₁ + I reconstruction the first component came out as 0
    # or a ~1e-16 roundoff floor; the chain exponential gives e^{λ} to
    # relative accuracy. Elementwise comparison on purpose.
    λ = [-240.0, -1.0]
    problem = IVP(LRHS(Diagonal(λ)), [1.0, 1.0], (0.0, 1.0))
    solution = solve(problem, ETD4RK(h = 0.25))
    @test all(isapprox.(solution.u[end], exp.(λ); rtol = 1e-10))
end

# A stable symmetric operator and a manufactured semilinear problem, shared
# by the testsets below. The problem u' = λu + f(u, t) with
# f(u, t) = -sin(t) - λcos(t) + (u - cos(t))² has exact solution u = cos(t)
# and a GENUINELY nonlinear (quadratic) non-stiff part, so all stages and
# couplings are exercised.
exprk_L = [-2.0 1.0 0.0; 1.0 -3.0 1.0; 0.0 1.0 -2.5]
exprk_u0 = [1.0, -0.5, 0.25]
exprk_solvers = (LawsonEuler, NorsettEuler, ETD2RK, ETD3RK, ETD4RK, Lawson4, Krogstad, HochbruckOstermann4)

@testset "exact linear propagation (V = eᶻ for every scheme)" begin
    problem = IVP(LRHS(exprk_L), exprk_u0, (0.0, 1.0))
    for maker in exprk_solvers
        solution = solve(problem, maker(h = 0.25))
        @test solution.u[end] ≈ exp(exprk_L) * exprk_u0  rtol = 1e-12
    end
end

@testset "NorsettEuler is exact for constant forcing at any h" begin
    g = [0.5, -1.0, 2.0]
    problem = IVP(LRHS(exprk_L, t -> g), exprk_u0, (0.0, 1.0))
    solution = solve(problem, NorsettEuler(h = 0.5))
    exact = exp(exprk_L) * exprk_u0 + exprk_L \ ((exp(exprk_L) - I) * g)
    @test solution.u[end] ≈ exact  rtol = 1e-11
end

@testset "classical convergence orders (manufactured nonlinear problem)" begin
    λ = -2.0
    f(u, t) = @. -sin(t) - λ * cos(t) + (u - cos(t))^2
    problem = IVP(SRHS(λ, f), [1.0], (0.0, 2.0))
    exact(t) = [cos(t)]
    for (maker, p) in ((LawsonEuler, 1), (NorsettEuler, 1), (ETD2RK, 2),
                       (ETD3RK, 3), (ETD4RK, 4), (Lawson4, 4),
                       (Krogstad, 4), (HochbruckOstermann4, 4))
        ps = sweep_params(p)
        order = convergence_order(h -> maker(h = h), problem, exact;
                                  h0 = ps.h0, refinements = ps.refinements)
        @test order > p - ps.below
    end
end

@testset "complex Diagonal (spectral) path" begin
    λs = ComplexF64[-1.0 + 2.0im, -3.0 - 1.0im, -0.5 + 0.0im]
    fc(u, t) = @. 0.1 * u^2
    u0c = ComplexF64[0.5, 0.2 - 0.1im, -0.3]
    problem = IVP(SRHS(Diagonal(λs), RHS(fc; iscomplex = true)), u0c, (0.0, 1.0))
    # Cross-family reference: classical RK4 at tiny h on the FULL right-hand side.
    full(u, t) = Diagonal(λs) * u .+ fc(u, t)
    reference = solve(IVP(RHS(full; iscomplex = true), u0c, (0.0, 1.0)), RK4(h = 1e-4))
    solution = solve(problem, ETD4RK(h = 1e-2))
    @test solution.u[end] ≈ reference.u[end]  rtol = 1e-7
end

@testset "sparse diagonal operators take the Diagonal fast path" begin
    λs = [-1.0, -2.0, -3.0]
    f2(u, t) = @. u^2 / 10
    u0 = [0.5, 0.3, 0.1]
    sol_sparse = solve(IVP(SRHS(spdiagm(0 => λs), f2), u0, (0.0, 1.0)), ETD4RK(h = 0.1))
    sol_diag   = solve(IVP(SRHS(Diagonal(λs),     f2), u0, (0.0, 1.0)), ETD4RK(h = 0.1))
    @test sol_sparse.u[end] ≈ sol_diag.u[end]  rtol = 1e-14
end

@testset "error paths" begin
    # A plain nonlinear problem has no L to exponentiate:
    nonlinear = IVP(RHS((u, t) -> @. -u^3), [1.0], (0.0, 1.0))
    @test_throws ArgumentError solve(nonlinear, ETD4RK(h = 0.1))
    # A split problem whose STIFF part is nonlinear is equally out of scope:
    # (plain broadcast dots here: a bare `@.` mid-argument-list is greedy and
    # would swallow the comma and the second function)
    splitnl = IVP(SRHS((u, t) -> -u .^ 3, (u, t) -> u ./ 10), [1.0], (0.0, 1.0))
    @test_throws ArgumentError solve(splitnl, ETD4RK(h = 0.1))
    # The φ-coefficients bind to z = h⋅L, so h must be strictly positive at
    # cache construction:
    okproblem = IVP(SRHS(-1.0, (u, t) -> @. u / 10), [1.0], (0.0, 1.0))
    @test_throws ArgumentError NSDEBase.initialize_cache(okproblem, ETD4RK())
    # ... and adaptivity is structurally rejected at solver construction:
    @test_throws ArgumentError EXPRK(ETD4RK(h = 0.1).tableau, 0.1, AP())
end

@testset "cache reuse (the Parareal path) is deterministic" begin
    λ = -2.0
    f(u, t) = @. -sin(t) - λ * cos(t) + (u - cos(t))^2
    problem = IVP(SRHS(λ, f), [1.0], (0.0, 2.0))
    solver = ETD4RK(h = 0.1)
    cache = NSDEBase.initialize_cache(problem, solver)
    solution = NSDEBase.initialize_solution(problem, solver)
    NSDEBase.solve!(cache, solution, problem, solver)
    first_run = copy(solution.u[end])
    NSDEBase.solve!(cache, solution, problem, solver) # reuse cache AND solution
    @test solution.u[end] == first_run
end
