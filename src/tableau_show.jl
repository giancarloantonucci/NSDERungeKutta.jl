# NSDERungeKutta/src/tableau_show.jl
#
# Butcher-array REPL printing — the survivor of the old to-do list's
# "Improve ButcherTableau printing" (its PrettyTables/AbstractTrees sibling
# is declined: this is 40 dependency-free lines). Coefficients stored as
# floats are shown as the exact rationals they are whenever a
# tolerance-guarded `rationalize` recovers one with a sane denominator
# (RK4's 1/6, DP54's 9017/3168); genuinely irrational entries (SDIRK's
# 1 ± √2/2 family) fall back to rounded floats rather than to the
# thousand-digit fractions naive rationalisation produces.
#
# Wire-up: `include("tableau_show.jl")` in src/NSDERungeKutta.jl after the
# tableau definitions. The MIME"text/plain" method takes REPL precedence
# over NSDEBase's generic AbstractObject show, which remains the compact
# fallback.
#
# Example:
#     julia> RK4(h = 1e-2).tableau
#     ButcherTableau (s = 4, p = 4):
#       0 │   0   0   0   0
#     1/2 │ 1/2   0   0   0
#     1/2 │   0 1/2   0   0
#       1 │   0   0   1   0
#     ────┼────────────────
#       4 │ 1/6 1/3 1/3 1/6

function _tableauentry(x::Real)
    isinteger(x) && return string(Int(x))
    r = rationalize(float(x); tol = √eps(float(x)))
    if abs(denominator(r)) ≤ 10_000 && isapprox(float(r), float(x); atol = 4eps(float(x)))
        return string(numerator(r), "/", denominator(r))
    end
    return string(round(float(x); sigdigits = 6))
end

function Base.show(io::IO, ::MIME"text/plain", tableau::ButcherTableau)
    s = tableau.s
    hasembedded = tableau.d !== nothing
    left = vcat(_tableauentry.(tableau.c), string(tableau.p),
                hasembedded ? [string(tableau.q)] : String[])
    rows = [[_tableauentry(tableau.A[i, j]) for j = 1:s] for i = 1:s]
    push!(rows, _tableauentry.(tableau.b))
    hasembedded && push!(rows, _tableauentry.(tableau.d))

    lw = maximum(length, left)
    cw = [maximum(length(rows[i][j]) for i in eachindex(rows)) for j = 1:s]

    print(io, "ButcherTableau (s = ", s, ", p = ", tableau.p)
    hasembedded && print(io, ", q = ", tableau.q)
    println(io, "):")
    rule = "─"^(lw + 1) * "┼" * "─"^(sum(cw) + s)
    for (i, row) in enumerate(rows)
        i == s + 1 && println(io, rule)          # c | A above, orders | weights below
        i == s + 2 && println(io, rule)          # embedded weights get their own rule
        print(io, lpad(left[i], lw), " │ ")
        println(io, join(lpad.(row, cw), " "))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Testset — paste into test/runtests.jl:
#
# @testset "ButcherTableau printing" begin
#     out = sprint(show, MIME("text/plain"), RK4(h = 1.0).tableau)
#     @test occursin("s = 4, p = 4", out)
#     @test occursin("1/6", out) && occursin("1/2", out)   # exact rationals recovered
#     @test occursin("┼", out)                             # the Butcher rule
#     @test length(findall("│", out)) == 5                 # 4 stage rows + weights
#
#     emb = sprint(show, MIME("text/plain"), DP54(h = 1.0).tableau)
#     @test occursin("q = ", emb)                          # embedded order in header
#     @test occursin("9017/3168", emb)                     # a DP54 coefficient, exactly
#     @test length(findall("┼", emb)) == 2                 # weights AND embedded rules
# end
