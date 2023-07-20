hairernorm(v) = sqrt(sum(abs2, v) / length(v))

# Kahan's compensated summation for improved numerical accuracy
function kahansum(a::AbstractFloat, b::AbstractFloat, e::Ref{<:AbstractFloat})
    corrected_b = b + e[]
    sum = a + corrected_b
    compensation = corrected_b - (sum - a)
    e[] = compensation
    return sum
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
