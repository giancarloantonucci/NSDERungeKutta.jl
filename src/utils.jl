Vector{u_T}(undef, n, d) where u_T = Vector{u_T}[Vector{u_T}(undef, d) for i = 1:n]
Vector{u_T}(undef, n₂, n₁, d) where u_T = Vector{Vector{u_T}}[Vector{u_T}(undef, n₁, d) for i = 1:n₂]

# zero!(v::AbstractVector{<:Number}) = fill!(v, zero(eltype(v)))
# function zero!(v::AbstractVector{<:AbstractVector{<:Number}})
#     for i in eachindex(v)
#         zero!(v[i])
#     end
#     return v
# end

function norm!(v::AbstractVector{<:AbstractVector{<:Number}})
    r = zero(eltype(v))
    for i in eachindex(v)
        r += norm(v[i])
    end
    return r
end

hairernorm(v) = sqrt(sum(v.^2) / length(v))

function +ₖ(y₀::AbstractFloat, tup::Tuple{AbstractFloat, Ref{<:AbstractFloat}})
    x₀, err = tup
    X₀ = x₀ + err[]
    y₁ = y₀ + X₀
    X̂₀ = y₁ - y₀
    err[] = X₀ - X̂₀
    return y₁
end
