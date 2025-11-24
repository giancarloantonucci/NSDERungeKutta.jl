# NSDERungeKutta/src/utils.jl

hairernorm(v) = sqrt(sum(abs2, v) / length(v))

"""
    compensated_sum(sum::T, addend::T, error::Ref{T}) where T<:AbstractFloat

Adds `addend` to `sum` using the Kahan-Babuška-Neumaier algorithm to minimize 
floating-point round-off error. The accumulated error is stored in the `error` Ref.
"""
@inline function compensated_sum(sum::T, addend::T, error::Ref{T}) where T<:AbstractFloat
    # 1. Recover the previous compensation
    v = addend + error[]
    
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
