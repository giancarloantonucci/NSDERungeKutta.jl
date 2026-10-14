# Butcher-array REPL display.

function _tableauentry(x::Real)
    isinteger(x) && return string(Int(x))
    r = rationalize(float(x); tol=(√eps(float(x))))
    if abs(denominator(r)) ≤ 10_000 && isapprox(float(r), float(x); atol=4eps(float(x)))
        return string(numerator(r), "/", denominator(r))
    end
    return string(round(float(x); sigdigits=6))
end

function Base.show(io::IO, ::MIME"text/plain", tableau::ButcherTableau)
    s = tableau.s
    hasembedded = tableau.d !== nothing
    left = vcat(
        _tableauentry.(tableau.c), string(tableau.p), hasembedded ? [string(tableau.q)] : String[]
    )
    rows = [[_tableauentry(tableau.A[i, j]) for j in 1:s] for i in 1:s]
    push!(rows, _tableauentry.(tableau.b))
    hasembedded && push!(rows, _tableauentry.(tableau.d))

    lw = maximum(length, left)
    cw = [maximum(length(rows[i][j]) for i in eachindex(rows)) for j in 1:s]

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
