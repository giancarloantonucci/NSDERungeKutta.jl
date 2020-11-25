@doc raw"""
    ℛ(z::Number, solver::RungeKuttaSolver)::Number

returns the stability function of a `RungeKuttaSolver`:
```math
    R(z) = \frac{\det(I - z(A - eb^\intercal))}{\det(I - zA)}.
```
"""
function ℛ(z::Number, tableau::ButcherTableau)
    @↓ A, b, s = tableau
    e = ones(s)
    return det(I - z * (A - e * b')) / det(I - z * A)
end
