@doc raw"""
    ButcherTableau{A_T, b_T, c_T, s_T, p_T, d_T, q_T}

defines a constructor for the Butcher tableau of a `RungeKuttaSolver`.

---

    ButcherTableau(A, b, c, s, p[, d, q])

returns a `ButcherTableau` with:
- `A :: AbstractMatrix` : matrix of coefficients.
- `b :: AbstractVector` : vector of weights.
- `c :: AbstractVector` : vector of nodes.
- `s :: Integer` : number of stages.
- `p :: Integer` : order of accuracy.
- `d :: AbstractVector` : embedding's vector of weights.
- `q :: Integer` : embedding's order of accuracy.

---

    ButcherTableau(tableau::AbstractMatrix)

returns a `ButcherTableau` from a matrix `tableau` with structure:
```math
\begin{array}{c|c}
    c & A \\
    \hline
    p & b \\
    q & d
\end{array}
```
"""
struct ButcherTableau{A_T, b_T, c_T, s_T, p_T, d_T, q_T}
    A::A_T
    b::b_T
    c::c_T
    s::s_T
    p::p_T
    d::d_T
    q::q_T
end

function ButcherTableau(A, b, c, s, p)
    d = nothing
    q = nothing
    return ButcherTableau(A, b, c, s, p, d, q)
end

function ButcherTableau(tableau)
    nrows, ncols = size(tableau)
    s = ncols - 1
    A = tableau[1:s, 2:ncols]
    b = tableau[s+1, 2:ncols]
    c = tableau[1:s, 1]
    p = convert(Integer, tableau[s+1, 1])
    if nrows == ncols
        return ButcherTableau(A, b, c, s, p)
    else
        d = tableau[nrows, 2:ncols]
        q = convert(Integer, tableau[nrows, 1])
        return ButcherTableau(A, b, c, s, p, d, q)
    end
end

Base.summary(io::IO, tableau::ButcherTableau) = print(io, "ButcherTableau")

function Base.show(io::IO, tableau::ButcherTableau)
    print(io, "ButcherTableau:")
    pad = get(io, :pad, "")
    names = propertynames(tableau)
    N = length(names)
    for (n, name) in enumerate(names)
        field = getproperty(tableau, name)
        if field !== nothing
            print(io, "\n", pad, "   ‣ " * string(name) * " ≔ ")
            show(IOContext(io, :pad => "   "), field)
        end
    end
end

"""
    StepSize{h_T}

returns a constructor containing the step-size of a `RungeKuttaSolver`.

---

    StepSize(; h)

returns a `StepSize` with:
- `h :: Real` : step-size.
"""
mutable struct StepSize{h_T}
    h::h_T
end

Base.summary(io::IO, stepsize::StepSize) = print(io, "StepSize")

function Base.show(io::IO, stepsize::StepSize)
    print(io, "StepSize:")
    pad = get(io, :pad, "")
    names = propertynames(stepsize)
    N = length(names)
    for (n, name) in enumerate(names)
        field = getproperty(stepsize, name)
        if field !== nothing
            print(io, "\n", pad, "   ‣ " * string(name) * " ≔ ")
            show(IOContext(io, :pad => "   "), field)
        end
    end
end

# Base.summary(io::IO, stepsize::StepSize) = print(io, stepsize.h)
#
# function Base.show(io::IO, stepsize::StepSize)
#     print(io, stepsize.h)
#     newline = get(io, :newline, "\n")
#     print(io, newline)
# end

"""
    AdaptiveParameters{δ_T, ϵ_T, K_T}

returns a constructor containing the parameters of an adaptive `RungeKuttaSolver`.

---

    AdaptiveParameters(; δ = 0.0, ϵ = 1e-5, K = 100)

returns an `AdaptiveParameters` with:
- `δ :: Real`    : absolute tolerance.
- `ϵ :: Real`    : relative tolerance.
- `K :: Integer` : maximum number of iterations.
"""
mutable struct AdaptiveParameters{δ_T, ϵ_T, K_T}
    δ::δ_T
    ϵ::ϵ_T
    K::K_T
end

AdaptiveParameters(; δ = 0.0, ϵ = 1e-5, K = 100) = AdaptiveParameters(δ, ϵ, K)

function adaptive_step!(solution, solver, cache)
    @↓ u = solution
    @↓ h = solver.stepsize
    @↓ n, m, v = cache
    if solver.adaptive isa AdaptiveParameters
        k = solution.k isa Nothing ? cache.k : solution.k[n]
        @↓ s, b, p, d, q = solver.tableau
        @↓ δ, ϵ, K = solver.adaptive
        zero!(v)
        for i = 1:s
            @. v += (b[i] - d[i]) * k[i]
        end
        @. v /= δ + max(abs(u[n]), abs(u[n+1])) * ϵ
        error = norm(v) / √length(v)
        power = - 1 / (min(p, q) + 1)
        factor = 0.9 * error ^ power
        h *= max(0.5, min((m == 1 ? 2.0 : 1.0), factor))
        if error < 1 || m ≥ K
            n += 1
        else
            m += 1
        end
    else
        n += 1
    end
    @↑ cache = n, m
    @↑ solver.stepsize = h
end

Base.summary(io::IO, adaptive::AdaptiveParameters) = print(io, "AdaptiveParameters")

function Base.show(io::IO, adaptive::AdaptiveParameters)
    print(io, "AdaptiveParameters:")
    pad = get(io, :pad, "")
    names = propertynames(adaptive)
    N = length(names)
    for (n, name) in enumerate(names)
        field = getproperty(adaptive, name)
        if field !== nothing
            print(io, "\n", pad, "   ‣ " * string(name) * " ≔ ")
            show(IOContext(io, :pad => "   "), field)
        end
    end
end

"""
    NewtonParameters{ϵ_T, K_T}

returns a constructor containing the parameters of simplified Newton in `ImplicitRungeKuttaSolver`.

---

    NewtonParameters(; ϵ = 1e-3, K = 10)

returns a `NewtonParameters` with:
- `ϵ :: Real`    : relative tolerance
- `K :: Integer` : maximum number of iterations.
"""
mutable struct NewtonParameters{ϵ_T, K_T}
    ϵ::ϵ_T
    K::K_T
end

NewtonParameters(; ϵ = 1e-3, K = 10) = NewtonParameters(ϵ, K)

Base.summary(io::IO, newton::NewtonParameters) = print(io, "NewtonParameters")

function Base.show(io::IO, newton::NewtonParameters)
    print(io, "NewtonParameters:")
    pad = get(io, :pad, "")
    names = propertynames(newton)
    N = length(names)
    for (n, name) in enumerate(names)
        field = getproperty(newton, name)
        if field !== nothing
            print(io, "\n", pad, "   ‣ " * string(name) * " ≔ ")
            show(IOContext(io, :pad => "   "), field)
        end
    end
end
