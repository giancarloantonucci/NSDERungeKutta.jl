@doc raw"""
    ButcherTableau <: AbstractButcherTableau

A composite type for the Butcher tableau of a Runge-Kutta solver:
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
ButcherTableau(tableau::AbstractMatrix{ℝ}) where ℝ<:Real
```

## Arguments
- `A :: AbstractMatrix{ℝ} where ℝ<:Real` : matrix of coefficients.
- `b :: AbstractVector{ℝ} where ℝ<:Real` : vector of weights.
- `c :: AbstractVector{ℝ} where ℝ<:Real` : vector of nodes.
- `s :: Integer` : number of stages.
- `p :: Integer` : order of accuracy.
- `d :: AbstractVector{ℝ} where ℝ<:Real` : embedding's vector of weights.
- `q :: Integer` : embedding's order of accuracy.

# Functions
[`butchertableau`](@ref) : return matrix of parameters.
"""
struct ButcherTableau{A_T<:(AbstractMatrix{ℝ} where ℝ<:Real), b_T<:(AbstractVector{ℝ} where ℝ<:Real), c_T<:(AbstractVector{ℝ} where ℝ<:Real), s_T<:Integer, p_T<:Integer, d_T<:(Union{AbstractVector{ℝ}, Nothing} where ℝ<:Real), q_T<:Union{Integer, Nothing}} <: AbstractButcherTableau
    A::A_T
    b::b_T
    c::c_T
    s::s_T
    p::p_T
    d::d_T
    q::q_T
end
ButcherTableau(A::AbstractMatrix{ℝ}, b::AbstractVector{ℝ}, c::AbstractVector{ℝ}, s::Integer, p::Integer) where ℝ<:Real = ButcherTableau(A, b, c, s, p, nothing, nothing)
function ButcherTableau(tableau::AbstractMatrix{ℝ}) where ℝ<:Real
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

#####
##### Functions
#####

"""
    butchertableau(solver::AbstractRungeKuttaSolver) :: AbstractMatrix{<:Real}

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
