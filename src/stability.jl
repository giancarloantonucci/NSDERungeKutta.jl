@doc raw"""
    ℛ(z::Number, tableau::AbstractButcherTableau) :: Number
    ℛ(z::Number, solver::AbstractRungeKuttaSolver) :: Number

returns the stability function of `solver`:
```math
    R(z) = \frac{\det(I - z(A - \mathbb{1}b^\intercal))}{\det(I - zA)}.
```
"""
function ℛ(z::Number, tableau::AbstractButcherTableau)
    @↓ A, b, s = tableau
    𝟙 = ones(s)
    return det(I - z * (A - 𝟙 * b')) / det(I - z * A)
end
ℛ(z::Number, solver::AbstractRungeKuttaSolver) = ℛ(z, solver.tableau)

"""
    ℛ(Z::AbstractMatrix, tableau::AbstractButcherTableau) :: AbstractMatrix
    ℛ(Z::AbstractMatrix, solver::AbstractRungeKuttaSolver) :: AbstractMatrix

returns the stability function of `solver`.
"""
function ℛ(Z::AbstractMatrix, tableau::AbstractButcherTableau)
    @↓ A, b, s = tableau
    𝟙 = ones(s)
    tmp = kron(𝟙, Z)
    tmp = (I - kron(A, Z)) \ tmp
    tmp = kron(b', Matrix(1.0I, size(Z)...)) * tmp
    return I + tmp
end
ℛ(Z::AbstractMatrix, solver::AbstractRungeKuttaSolver) = ℛ(Z, solver.tableau)
