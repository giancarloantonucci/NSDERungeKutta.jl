"""
    NewtonParameters <: AbstractNewtonParameters

A composite type for the parameters of simplified Newton.

# Constructors
```julia
NewtonParameters(; rtol=1e-3, nits=10)
```

## Arguments
- `rtol :: Real` : relative tolerance
- `nits :: Integer` : maximum number of iterations allowed
"""
struct NewtonParameters{rtol_T<:Real, nits_T<:Integer} <: AbstractNewtonParameters
    rtol::rtol_T
    nits::nits_T
end

NewtonParameters(; rtol::Real=1e-3, nits::Integer=10) = NewtonParameters(rtol, nits)
