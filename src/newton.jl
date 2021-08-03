"""
    NewtonParameters

A composite type for the parameters of the parameters of simplified Newton in [`ImplicitRungeKuttaSolver`](@ref).

# Constructors
```julia
NewtonParameters(; ϵ = 1e-3, K = 10)
```

# Arguments
- `ϵ :: Real`    : relative tolerance
- `K :: Integer` : maximum number of iterations.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
mutable struct NewtonParameters{ϵ_T, K_T}
    ϵ::ϵ_T
    K::K_T
end

NewtonParameters(; ϵ = 1e-3, K = 10) = NewtonParameters(ϵ, K)

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

"""
    show(io::IO, newton::NewtonParameters)

prints a full description of `newton` and its contents to a stream `io`.
"""
Base.show(io::IO, newton::NewtonParameters) = NSDEBase._show(io, newton)

"""
    summary(io::IO, newton::NewtonParameters)

prints a brief description of `newton` to a stream `io`.
"""
Base.summary(io::IO, newton::NewtonParameters) = NSDEBase._summary(io, newton)
