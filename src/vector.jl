Vector{u_T}(undef, n, d) where u_T = Vector{u_T}[Vector{u_T}(undef, d) for i = 1:n]
Vector{u_T}(undef, n₂, n₁, d) where u_T = Vector{Vector{u_T}}[Vector{u_T}(undef, n₁, d) for i = 1:n₂]

zero!(v::AbstractVector) = fill!(v, zero(eltype(v)))
function zero!(v::AbstractVector{<:AbstractVector})
    for i in eachindex(v)
        zero!(v[i])
    end
    return v
end

function norm!(v::AbstractVector{<:AbstractVector})
    r = zero(eltype(v))
    for i in eachindex(v)
        r += norm(v[i])
    end
    return r
end
