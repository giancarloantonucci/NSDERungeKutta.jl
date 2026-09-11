# NSDERungeKutta/src/exprk/tableau.jl

@doc raw"""
    ExponentialTableau <: AbstractRungeKuttaParameters

A composite type for the coefficient functions of a one-step exponential
Runge-Kutta method for ``u' = Lu + g(t) + f_\text{ns}(u, t)``. Unlike a
[`ButcherTableau`](@ref), whose entries are numbers, an exponential tableau's
entries are OPERATORS — linear combinations of ``\varphi``-functions and
exponentials evaluated at ``z = hL`` — so the tableau stores a coefficient
FUNCTION `κ` and the solver cache evaluates it once per solve.

Following EXPINT's scheme-file convention (with `nothing` in place of Matlab's
`[]`), `κ(z)` must return the tuple `(U, V, A, B)` where, for `s` stages,
- `U :: Vector` (length `s`) : the operators applied to `u[n]` in each stage,
- `V` : the operator applied to `u[n]` in the update (``e^z`` for every
  scheme implemented here),
- `A :: Matrix` (`s × s`, strictly lower triangular, `nothing` for absent
  entries) : the stage-coupling operators,
- `B :: Vector` (length `s`, `nothing` for absent entries) : the weights,

so that a step reads ``U_i = U_i u_n + h \sum_{j<i} A_{ij} k_j``,
``k_i = f_\text{ns}(U_i, t_n + c_i h) + g(t_n + c_i h)`` and
``u_{n+1} = V u_n + h \sum_i B_i k_i``.

# Arguments
- `name :: Symbol` : the method's name.
- `p :: Integer` : classical (non-stiff) order.
- `q :: Integer` : stiff order in the sense of Hochbruck–Ostermann, as stated
  in EXPINT's scheme headers. For semilinear stiff PDEs the observed order is
  governed by `q`, not `p`.
- `s :: Integer` : number of stages.
- `c :: AbstractVector{<:Real}` : quadrature nodes.
- `κ :: Function` : the coefficient function described above.
"""
struct ExponentialTableau{c_T<:AbstractVector{<:Real}, κ_T<:Function} <: AbstractRungeKuttaParameters
    name :: Symbol
    p :: Int
    q :: Int
    s :: Int
    c :: c_T
    κ :: κ_T
    function ExponentialTableau(name::Symbol, p::Integer, q::Integer, s::Integer, c::c_T, κ::κ_T) where {c_T<:AbstractVector{<:Real}, κ_T<:Function}
        length(c) == s || throw(ArgumentError("`ExponentialTableau` needs `length(c) == s`."))
        return new{c_T, κ_T}(name, Int(p), Int(q), Int(s), c, κ)
    end
end
