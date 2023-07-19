hairernorm(v) = sqrt(sum(v.^2) / length(v))

# Kahan's (compensated) summation
function kahan_sum(a::AbstractFloat, b::AbstractFloat, e::Ref{AbstractFloat})
    b1 = b + e[]
    c = a + b1
    b2 = c - a
    e[] = b1 - b2
    return c
end

# interpolation methods
function linearspline(x, xᵢ₋₁, xᵢ, yᵢ₋₁, yᵢ)
    hᵢ = xᵢ - xᵢ₋₁
    aᵢ₋₁ = (xᵢ - x) / hᵢ
    aᵢ = (x - xᵢ₋₁) / hᵢ
    y = @. aᵢ₋₁ * yᵢ₋₁ + aᵢ * yᵢ
    return y
end

function hermitecubicspline(x, xᵢ₋₁, xᵢ, yᵢ₋₁, yᵢ, dyᵢ₋₁, dyᵢ)
    hᵢ = xᵢ - xᵢ₋₁
    c₀ = yᵢ₋₁
    c = dyᵢ₋₁
    c₂ = @. (3 * (yᵢ - yᵢ₋₁) / hᵢ - (dyᵢ + 2 * dyᵢ₋₁)) / hᵢ
    c₃ = @. ((dyᵢ + dyᵢ₋₁) - 2 * (yᵢ - yᵢ₋₁) / hᵢ) / hᵢ^2
    y = @. c₀ + c * (x - xᵢ₋₁) + c₂ * (x - xᵢ₋₁)^2 + c₃ * (x - xᵢ₋₁)^3
    return y
end
