@doc raw"""
    ButcherTableau <: AbstractButcherTableau

A composite type for the Butcher tableau of an [`AbstractRungeKuttaSolver`](@ref):
```math
\begin{array}{c|c}
    c & A \\
    \hline
    p & b \\
    q & d
\end{array}
```

# Constructors
```julia
ButcherTableau(A, b, c, s, p[, d, q])
ButcherTableau(tableau)
```

# Arguments
- `A :: AbstractMatrix` : matrix of coefficients.
- `b :: AbstractVector` : vector of weights.
- `c :: AbstractVector` : vector of nodes.
- `s :: Integer` : number of stages.
- `p :: Integer` : order of accuracy.
- `d :: AbstractVector` : embedding's vector of weights.
- `q :: Integer` : embedding's order of accuracy.
- `tableau :: AbstractMatrix` : from which all other fields can be extracted.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct ButcherTableau{A_T, b_T, c_T, s_T, p_T, d_T, q_T} <: AbstractButcherTableau
    A::A_T
    b::b_T
    c::c_T
    s::s_T
    p::p_T
    d::d_T
    q::q_T
end

function ButcherTableau(A::AbstractMatrix, b::AbstractVector, c::AbstractVector, s::Integer, p::Integer)
    d = nothing
    q = nothing
    return ButcherTableau(A, b, c, s, p, d, q)
end

function ButcherTableau(tableau::AbstractMatrix)
    nrows, ncols = size(tableau)
    s = ncols - 1
    A = tableau[1:s, 2:ncols]
    b = tableau[s+1, 2:ncols]
    c = tableau[1:s, 1]
    p = convert(Integer, tableau[s+1, 1])
    if nrows == ncols
        return ButcherTableau(A, b, c, s, p)
    else
        d = tableau[nrows, 2:ncols]
        q = convert(Integer, tableau[nrows, 1])
        return ButcherTableau(A, b, c, s, p, d, q)
    end
end
