# NSDERungeKutta/src/stability.jl

@doc raw"""
    stability_function(z::Number, tableau::AbstractButcherTableau) :: Number
    stability_function(z::Number, solver::AbstractRungeKuttaSolver) :: Number

Computes the scalar stability function $R(z) = 1 + z b^\top (I - zA)^{-1} \mathbb{1}$.
"""
function stability_function(z::Number, tableau::AbstractButcherTableau)
    @↓ A, b, s = tableau

    if iszero(z)
        return one(z)
    end

    # Compute w = (I - zA)⁻¹ * 1
    if is_strictly_lower_triangular(A)
        # Forward substitution for unit-lower (I - zA)
        # Promote type to ensure we don't put Floats into an Int array
        T = promote_type(typeof(z), eltype(A))
        w = ones(T, s)

        for i in 1:s
            sum_Aj = zero(T)
            @inbounds for j in 1:i-1
                sum_Aj += A[i, j] * w[j]
            end
            w[i] = 1 + z * sum_Aj
        end
    else
        I_minus_zA = I - z * A
        w = I_minus_zA \ ones(typeof(z), s)
    end

    # R(z) = 1 + z * bᵀ w
    return one(z) + z * dot(b, w)
end

stability_function(z::Number, solver::AbstractRungeKuttaSolver) = stability_function(z, solver.tableau)

@doc raw"""
    stability_function(Z::AbstractMatrix, tableau::AbstractButcherTableau) :: AbstractMatrix
    stability_function(Z::AbstractMatrix, solver::AbstractRungeKuttaSolver) :: AbstractMatrix

Computes the matrix stability function $R(Z) = I + (b^\top \otimes Z) (I \otimes I - A \otimes Z)^{-1} (\mathbb{1} \otimes I)$.
This builds and solves a Kronecker system of size $sN \times sN$; avoid for large systems.
"""
function stability_function(Z::AbstractMatrix, tableau::AbstractButcherTableau)
    @↓ A, b, s = tableau
    d = size(Z, 1)

    I_d = Matrix{eltype(Z)}(I, d, d)
    ones_s = ones(s)

    # Solve (I - A ⊗ Z) * Y = (1 ⊗ Z)
    LHS = I - kron(A, Z)
    RHS = kron(ones_s, Z)
    Y = LHS \ RHS
    
    # R(Z) = I + (bᵀ ⊗ I) Y
    return I_d + kron(b', I_d) * Y
end

stability_function(Z::AbstractMatrix, solver::AbstractRungeKuttaSolver) = stability_function(Z, solver.tableau)

"""
    is_strictly_lower_triangular(A::AbstractMatrix) :: Bool

Returns `true` if all entries on and above the main diagonal of `A` are zero.
"""
function is_strictly_lower_triangular(A::AbstractMatrix)
    for j in axes(A, 2), i in 1:j
        if !iszero(A[i,j])
            return false
        end
    end
    return true
end
