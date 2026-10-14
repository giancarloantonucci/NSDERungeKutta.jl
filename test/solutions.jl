@testset "dense output and interpolation" begin
    @testset "stage history aligns with accepted steps" begin
        problem = Logistic(0.3, (0.0, 1.0))
        solver = BogackiShampine(h=0.05, εᵣ=1e-6)
        solution = solve(problem, solver; dense=true)
        @test length(solution.k) == numtimesteps(solution) - 1 # the drift fix
        for tₚ in (0.13, 0.5, 0.87)
            @test norm(solution(tₚ, solver.tableau) - logistic_exact(tₚ)) < 1e-4
        end
    end

    @testset "spline fallbacks" begin
        solution = solve(dahlquist, RungeKutta4(h=0.1))
        exact = dahlquist_exact(0.55)
        err_linear = norm(solution(0.55) - exact)
        err_hermite = norm(solution(0.55, (u, t) -> u) - exact) # u′ = u
        @test err_linear < 5e-3
        @test err_hermite < 1e-5
        @test err_hermite < err_linear
        @test solution(-10.0) == solution.u[1]   # clamped below
        @test solution(99.0) == solution.u[end]  # clamped above
        @test solution(solution.t[3]) == solution.u[3] # exact at grid nodes (pins the segment lookup)
    end
end

@testset "IMEX rejects unsupported dense history before stepping" begin
    for make_solver in ierk_zoo
        solver = make_solver(h=0.1)
        @test_throws ArgumentError RungeKuttaSolution(split_lin, solver; dense=true)
        @test_throws ArgumentError solve(split_lin, solver; dense=true)
        # The supported path still integrates and interpolates normally.
        sol = solve(split_lin, solver)
        @test sol.k === nothing
        @test all(isfinite, sol(0.15))
    end
    # A caller may bypass construction by reusing a dense solution. Reject
    # that path too, before clearing its history or changing the cache.
    solver = IMEXEuler(h=0.1)
    cache = NSDERungeKutta.RungeKuttaCache(split_lin, solver)
    sol = solve(dahlquist, Euler(h=0.1); dense=true)
    before = deepcopy(sol)
    @test_throws ArgumentError NSDEBase.solve!(cache, sol, split_lin, solver)
    @test sol.u == before.u && sol.t == before.t && sol.k == before.k
    @test cache.n == 1
end

@testset "solution container" begin
    solution = solve(dahlquist, RungeKutta4(h=0.25))
    @test numvariables(solution) == 1
    @test numtimesteps(solution) == length(solution.t)
    @test length(extract(solution, 1)) == numtimesteps(solution)

    # setindex! with a solution slice (the fixed `@↓ … ← …` path)
    solution[2] = solution[1]
    @test solution.u[2] == solution.u[1]
    @test solution.t[2] == solution.t[1]
    # tuple form
    solution[3] = ([9.9], 9.9)
    @test solution.u[3] == [9.9] && solution.t[3] == 9.9
end

@testset "dense slicing: stages live on intervals, not nodes" begin
    # A dense solution with N nodes stores N − 1 stage sets. The indexers used
    # to read `k[i]` at node i, so `sol[end]` threw a BoundsError.
    dense = solve(dahlquist, Euler(h=0.5); dense=true)   # t = 0, 0.5, 1
    N = length(dense)
    @test length(dense.k) == N - 1
    last = dense[end]                                          # used to throw
    @test last.u == [dense.u[end]] && last.t == [dense.t[end]] && last.k === nothing
    @test dense[1].k === nothing                               # a node owns no interval
    slice = dense[1:2]                                         # two nodes, one interval
    @test slice.u == dense.u[1:2] && length(slice.k) == 1 && slice.k[1] == dense.k[1]
    whole = dense[1:N]
    @test length(whole.k) == N - 1                             # `end` inside a range is fine too
    @test_throws ArgumentError dense[[1, N]]                   # non-adjacent: no stages between them
    sparse = solve(dahlquist, RungeKutta4(h=0.5))
    @test sparse[[1, N]].k === nothing                         # non-dense: any selection is fine
    # Writes into a DENSE solution: partial writes are refused BEFORE any
    # mutation (the adjoining intervals' stages could no longer match), a
    # whole-solution write with matching stages goes through.
    before = (copy.(dense.u), copy(dense.t), deepcopy(dense.k))
    @test_throws ArgumentError dense[2] = dense[1]
    @test_throws ArgumentError dense[2] = ([9.9], 9.9)
    @test_throws ArgumentError dense[1:2] = slice
    @test_throws ArgumentError dense[1:N] = sparse                # nodes only: no stages to go with them
    @test_throws DimensionMismatch dense[1:2] = dense[1:N]        # sizes disagree
    @test dense.u == before[1] && dense.t == before[2] && dense.k == before[3] # untouched by every refusal
    other = solve(Dahlquist(0.5, (0.0, 1.0); λ=-2.0), Euler(h=0.5); dense=true)
    dense[1:N] = other
    @test dense.u == other.u && dense.t == other.t && dense.k == other.k
    # Non-dense targets accept node writes as before:
    sparse[2] = sparse[1]
    @test sparse.u[2] == sparse.u[1]
    u1, u3 = copy(sparse.u[1]), copy(sparse.u[3])
    sparse[[1, 3]] = sparse[[3, 1]]                               # non-contiguous is fine without stages
    @test sparse.u[1] == u3 && sparse.u[3] == u1
end

@testset "RecipesBase recipes (headless)" begin
    # `apply_recipe` exercises the Plots-side recipe surface with no Plots
    # dependency — previously the one exported surface the suite never
    sol = solve(Lorenz([2.0, 3.0, -14.0], (0.0, 1.0)), RK4(h=1e-2))
    for attrs in (
        Dict{Symbol,Any}(),
        Dict{Symbol,Any}(:variables => 2),
        Dict{Symbol,Any}(:iscomplex => false, :skip => 10),
    )
        @test RecipesBase.apply_recipe(attrs, sol) isa Vector{RecipesBase.RecipeData}
    end
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(), NSDEBase._PhasePlot(sol)) isa
        Vector{RecipesBase.RecipeData}
    zsol = solve(Dahlquist([1.0 + 0.0im], (0.0, 1.0); λ=-0.5 + 1.0im), RK4(h=1e-2))
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:iscomplex => true), zsol) isa
        Vector{RecipesBase.RecipeData}
    # The step-size record recipes, including the trailing-bin guard:
    recording = DormandPrince54(h=0.1, εᵣ=1e-4, save_stepsizes=true)
    rsol = solve(Logistic(0.3, (0.0, 1.0)), recording)
    hs = rsol.stepsizes
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(), hs) isa Vector{RecipesBase.RecipeData}
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(), rsol.t[1:(end - 1)], hs) isa
        Vector{RecipesBase.RecipeData}
end
