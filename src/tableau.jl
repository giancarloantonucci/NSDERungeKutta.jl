@doc raw"""
    ButcherTableau <: AbstractButcherTableau

A composite type for the Butcher tableau of an [`AbstractRungeKuttaSolver`](@ref):
```math
\begin{array}{c|c}
    c & A \\
    \hline
    p & b^\intercal \\
    q & d^\intercal
\end{array}
```

# Constructors
```julia
ButcherTableau(A, b, c, s, p[, d, q])
ButcherTableau(tableau)
```

## Arguments
- `A :: AbstractMatrix{<:Real}` : matrix of coefficients.
- `b :: AbstractVector{<:Real}` : vector of weights.
- `c :: AbstractVector{<:Real}` : vector of nodes.
- `s :: Integer` : number of stages.
- `p :: Integer` : order of accuracy.
- `d :: AbstractVector{<:Real}` : embedding's vector of weights.
- `q :: Integer` : embedding's order of accuracy.
- `tableau :: AbstractMatrix` : matrix of parameters as in the definition.

# Functions
[`butchertableau`](@ref) : returns Butcher tableau as matrix.
"""
struct ButcherTableau{A_T<:AbstractMatrix{<:Real}, b_T<:AbstractVector{<:Real}, c_T<:AbstractVector{<:Real}, s_T<:Integer, p_T<:Integer, d_T<:Union{AbstractVector{<:Real}, Nothing}, q_T<:Union{Integer, Nothing}} <: AbstractButcherTableau
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

"""
    butchertableau(solver::AbstractRungeKuttaSolver)

returns the Butcher tableau of a `solver` as a matrix.
"""
function butchertableau(solver::AbstractRungeKuttaSolver)
    @↓ A, b, c, s, p, d, q = solver.tableau
    if d isa Nothing && q isa Nothing
        return [c A; p transpose(b)]
    else
        return [c A; p transpose(b); q transpose(d)]
    end
end
