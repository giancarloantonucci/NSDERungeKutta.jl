hairernorm(v) = sqrt(sum(v.^2) / length(v))

# Kahan's (compensated) summation
function +ₖ(a::AbstractFloat, tpl::Tuple{AbstractFloat, Ref{<:AbstractFloat}})
    b, ε = tpl
    b₁ = b + ε[]
    c = a + b₁
    b₂ = c - a
    ε[] = b₁ - b₂
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
