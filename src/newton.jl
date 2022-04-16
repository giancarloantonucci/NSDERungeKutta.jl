"""
    SimplifiedNewtonParameters <: AbstractSimplifiedNewtonParameters

A composite type for the parameters of simplified Newton.
    
# Constructors
```julia
SimplifiedNewtonParameters(; rtol=1e-3, nits=10)
```

## Arguments
- `rtol :: Real` : relative tolerance
- `nits :: Integer` : maximum number of iterations allowed
"""
struct SimplifiedNewtonParameters{rtol_T<:Real, nits_T<:Integer} <: AbstractSimplifiedNewtonParameters
    rtol::rtol_T
    nits::nits_T
end
SimplifiedNewtonParameters(; rtol::Real=1e-3, nits::Integer=10) = SimplifiedNewtonParameters(rtol, nits)
