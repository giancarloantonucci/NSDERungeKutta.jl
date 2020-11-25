"""
    ButcherTableau(A, b, c, s[, b̂]) -> ButcherTableau
    ButcherTableau(tableau) -> ButcherTableau

returns a constructor for the `ButcherTableau` of a `RungeKuttaSolver` with
- `A` : matrix of coefficients.
- `b` : vector of weights.
- `c` : vector of nodes.
- `s` : number of stages.
- `b̂` : vector of weights (for adaptive methods).
"""
struct ButcherTableau{A_T, b_T, c_T, s_T, b̂_T}
    A::A_T
    b::b_T
    c::c_T
    s::s_T
    b̂::b̂_T
end

ButcherTableau(A, b, c, s) = ButcherTableau(A, b, c, s, similar(b))

function ButcherTableau(tableau)
    nrows, ncols = size(tableau)
    s = ncols - 1
    A = tableau[1:s, 2:ncols]
    b = tableau[s+1, 2:ncols]
    c = tableau[1:s, 1]
    if nrows == ncols
        return ButcherTableau(A, b, c, s)
    else
        b̂ = tableau[nrows, 2:ncols]
        return ButcherTableau(A, b, c, s, b̂)
    end
end
