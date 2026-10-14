using RootedTrees

# Test all trees, not only the linear chains seen by Dahlquist. Work with
# unnormalised residuals so large symmetry factors do not hide defects.
@testset "Rooted-tree order conditions" begin
    for F in [erk_fixed; erk_embedded; dirk_zoo; irk_zoo]
        tab = F(h=0.1).tableau
        @test vec(sum(tab.A; dims=2)) ≈ tab.c atol=2e-13 rtol=2e-13
        for (name, b, p) in ((:primary, tab.b, tab.p), (:embedded, tab.d, tab.q))
            b === nothing && continue
            rk = RungeKuttaMethod(tab.A, b, tab.c)
            @testset "$F $name order $p" begin
                for q in 1:p, tree in RootedTreeIterator(q)
                    @test abs(elementary_weight(tree, rk) - 1/density(tree)) < 1e-12
                end
                @test any(
                    abs(elementary_weight(t, rk)-1/density(t)) > 1e-12 for
                    t in RootedTreeIterator(p+1)
                )
            end
        end
        if tab.b_dense !== nothing
            @test vec(sum(tab.b_dense; dims=1)) ≈ tab.b atol=2e-14
        end
    end
end

@testset "Continuous extensions satisfy tree conditions inside a step" begin
    for (F, p) in ((Euler, 1), (BogackiShampine, 3), (DormandPrince54, 4))
        tab = F(h=0.1).tableau
        for θ in (0.0, 0.13, 0.5, 0.87, 1.0)
            bθ = vec(sum(tab.b_dense .* [θ^j for j in axes(tab.b_dense, 1)]; dims=1))
            rk = RungeKuttaMethod(tab.A, bθ, tab.c)
            for q in 1:p, tree in RootedTreeIterator(q)
                @test abs(elementary_weight(tree, rk) - θ^q/density(tree)) < 2e-13
            end
        end
    end
    # Ralston's tableau rounded to eight digits passes a coarse convergence fit
    # but fails even second-order trees at the tolerance used above.
    A = [0.0 0 0 0; 0.4 0 0 0; 0.29697761 0.15875964 0 0; 0.21810040 -3.05096516 3.83286476 0]
    b = [0.17476028, -0.55148066, 1.20553560, 0.17118478]
    rounded = RungeKuttaMethod(A, b)
    @test any(
        abs(elementary_weight(t, rounded) - 1 / density(t)) > 1e-12 for t in RootedTreeIterator(2)
    )
end
