zero!(v::AbstractVector{<:Number}) = fill!(v, zero(eltype(v)))
function zero!(v::AbstractVector{<:AbstractVector{<:Number}})
    for i in eachindex(v)
        zero!(v[i])
    end
    return v
end

function norm!(v::AbstractVector{<:AbstractVector{<:Number}})
    r = zero(eltype(v))
    for i in eachindex(v)
        r += norm(v[i])
    end
    return r
end
