# NSDERungeKutta/src/utils.jl

"""
    compensated_sum(sum::T, addend::Real, error::Ref{T}) where T<:AbstractFloat

Adds `addend` to `sum` using the Kahan-Babuška-Neumaier algorithm to minimize 
floating-point round-off error. The step is converted to the problem time type `T`.
The accumulated error is stored in the `error` Ref.
"""
@inline function compensated_sum(sum::T, addend::Real, error::Ref{T}) where {T<:AbstractFloat}
    # 1. Recover the previous compensation
    v = T(addend) + error[]

    # 2. Perform the addition
    new_sum = sum + v

    # 3. Update the error (Neumaier)
    if abs(sum) >= abs(v)
        error[] = (sum - new_sum) + v
    else
        error[] = (v - new_sum) + sum
    end

    return new_sum
end

# Linear spline interpolation between two points
function linearspline(x, x_prev, x_curr, y_prev, y_curr)
    h = x_curr - x_prev
    a_prev = (x_curr - x) / h
    a_curr = (x - x_prev) / h
    y = @. a_prev * y_prev + a_curr * y_curr
    return y
end

# Cubic Hermite spline interpolation between two points
function hermitecubicspline(x, x_prev, x_curr, y_prev, y_curr, dy_prev, dy_curr)
    h = x_curr - x_prev
    c0 = y_prev
    c1 = dy_prev
    c2 = @. (3 * (y_curr - y_prev) / h - (dy_curr + 2 * dy_prev)) / h
    c3 = @. ((dy_curr + dy_prev) - 2 * (y_curr - y_prev) / h) / h^2
    y = @. c0 + c1 * (x - x_prev) + c2 * (x - x_prev)^2 + c3 * (x - x_prev)^3
    return y
end

"""
    directldiv!(M, v)

In-place left-division `v = M \\ v` for a factorisation `M`. Identical to
`LinearAlgebra.ldiv!(M, v)` except for factorisation types that ship without a
two-argument in-place method — notably `SparseArrays.CHOLMOD.Factor`, which
`factorize` returns for symmetric positive-definite SPARSE matrices (i.e.
every method-of-lines Laplacian fed to the direct-linear DIRK/IERK/IRK
branches as `I - hA⋅L`) and which supports only `\\` (checked through Julia 1.13). The
fallback allocates one vector per call; the per-step `factorize` sitting next
to every call site already allocates strictly more, so the hot-path cost is
unchanged in order. CHOLMOD's `\\` also handles a complex right-hand side
against a real factor (real/imaginary split), covering the `iscomplex`
spectral path.
"""
directldiv!(M, v) = ldiv!(M, v)

# CHOLMOD factors provide an allocating solve rather than two-argument ldiv!.
directldiv!(M::SparseArrays.CHOLMOD.Factor, v) = copyto!(v, M \ v)
