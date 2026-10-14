# NSDERungeKutta/src/exprk/cache.jl

mutable struct ExponentialRungeKuttaCache{
    n_T<:Integer,
    e_T<:Ref{<:AbstractFloat},
    v_T<:AbstractVector{<:Number},
    k_T<:AbstractVector{v_T},
    op_T,
    A_T<:AbstractMatrix,
    B_T<:AbstractVector,
    h_T<:Real,
} <: AbstractRungeKuttaCache
    n::n_T          # step counter
    m::n_T          # adaptive step counter (always advances by 1; kept for the shared solve! loop)
    e::e_T          # compensated summation error
    v::v_T          # stage-state accumulator / update accumulator
    w::v_T          # scratch for the forcing g(t)
    k::k_T          # stage values fₙₛ(Uᵢ, tᵢ) + g(tᵢ) at step `n` (UNSCALED: h is folded into the mul!s)
    # Coefficient operators bind this cache to L and h.
    U::Vector{op_T}
    V::op_T
    A::A_T
    B::B_T
    h::h_T         # step used to build the coefficient operators
    tail::Union{Nothing,ExponentialRungeKuttaCache}
end

# The stiff part must be explicitly linear: exponential integrators propagate
# e^{hL} exactly and quadrature the rest. Anything else fails loudly here, at
# cache construction, before a single step is taken.
function linearoperator(
    rhs::SplitRightHandSide{𝐿,𝑁}
) where {𝐿<:LinearRightHandSide,𝑁<:NonlinearRightHandSide}
    return rhs.fₛ.L
end
linearoperator(rhs::LinearRightHandSide) = rhs.L
function linearoperator(rhs::AbstractRightHandSide)
    return throw(
        ArgumentError(
            "Exponential Runge-Kutta solvers integrate u' = Lu + g(t) + fₙₛ(u, t) and need " *
            "the stiff part supplied as a `LinearRightHandSide`, i.e. `SRHS(L, fₙₛ)`, " *
            "`SRHS(LRHS(L, g), fₙₛ)`, or a plain `LRHS(L[, g])`. Got a `$(typeof(rhs))`.",
        ),
    )
end

# Sparse diagonal operators (e.g. spectral discretisations built with
# `spdiagm`) take the elementwise Diagonal fast path; general sparse operators
# are densified, since their φ-functions are dense anyway (see
# [`phifunctions`](@ref)).
operatorform(L::AbstractMatrix) = L
function operatorform(L::SparseArrays.AbstractSparseMatrixCSC)
    return isdiag(L) ? Diagonal(Vector(diag(L))) : Matrix(L)
end

function RungeKuttaCache(problem::AbstractInitialValueProblem, solver::ExponentialRungeKuttaSolver)
    @↓ u0, tspan = problem
    t0, tN = tspan
    @↓ tableau, stepsize = solver
    @↓ s, κ = tableau
    @↓ h = stepsize
    h > 0 || throw(
        ArgumentError(
            "`ExponentialRungeKuttaSolver` needs `h > 0` at cache construction: the φ-coefficients are precomputed at z = h⋅L.",
        ),
    )
    L = operatorform(linearoperator(problem.rhs))
    opU, opV, opA, opB = κ(h * L)
    Op = typeof(opV)
    n = m = 1
    e = Ref(zero(eltype(t0)))
    v = similar(u0)
    w = similar(u0)
    k = [similar(u0) for i in 1:s]
    return ExponentialRungeKuttaCache(
        n,
        m,
        e,
        v,
        w,
        k,
        Vector{Op}(opU),
        opV,
        Matrix{Union{Op,Nothing}}(opA),
        Vector{Union{Op,Nothing}}(opB),
        h,
        nothing,
    )
end
