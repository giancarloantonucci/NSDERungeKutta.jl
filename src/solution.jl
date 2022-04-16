"""
    RungeKuttaSolution <: AbstractRungeKuttaSolution

A composite type for an [`AbstractRungeKuttaSolution`](@ref) obtained using an [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
RungeKuttaSolution(u, t)
RungeKuttaSolution(problem, solver)
```

## Arguments
- `u :: AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number` : numerical solution.
- `t :: AbstractVector{ℝ} where ℝ<:Real` : time grid.

# Functions
- [`extract`](@ref) : extract variable.
- [`firstindex`](@ref) : first index.
- [`getindex`](@ref) : get value(s) and time.
- [`lastindex`](@ref) : last index.
- [`length`](@ref) : number of time steps.
- [`setindex!`](@ref) : set value(s) and time.
- [`size`](@ref) : number of variables and time steps.
"""
struct RungeKuttaSolution{u_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), t_T<:(AbstractVector{ℝ} where ℝ<:Real)} <: AbstractRungeKuttaSolution
    u::u_T
    t::t_T
end
function RungeKuttaSolution(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ h = solver.stepsize
    N = round(Integer, (tN - t0) / h) + 1
    # initialise u
    u0_T = eltype(u0)
    d = length(u0)
    u = Vector{u0_T}(undef, N, d)
    u[1] = u0
    # initialise t
    t0_T = eltype(t0)
    t = Vector{t0_T}(undef, N)
    t[1] = t0
    return RungeKuttaSolution(u, t)
end

#####
##### Methods
#####

function findindex(tₚ, t, u)
    N = length(t)
    if tₚ < t[1]
        return 2
    end
    for n in 2:N
        if t[n-1] ≤ tₚ ≤ t[n]
            return n
        end
    end
    if tₚ > t[N]
        return N
    end
end

function nearestneighbour(x, xᵢ₋₁, xᵢ, fᵢ₋₁, fᵢ)
    s = abs(xᵢ - x) > abs(x - xᵢ₋₁) ? fᵢ₋₁ : fᵢ
    return s
end

function linearspline(x, xᵢ₋₁, xᵢ, fᵢ₋₁, fᵢ)
    hᵢ = xᵢ - xᵢ₋₁
    aᵢ₋₁ = (xᵢ - x) / hᵢ
    aᵢ = (x - xᵢ₋₁) / hᵢ
    s = @. aᵢ₋₁ * fᵢ₋₁ + aᵢ * fᵢ
    return s
end

function hermitecubicspline(x, xᵢ₋₁, xᵢ, fᵢ₋₁, fᵢ, dfᵢ₋₁, dfᵢ)
    hᵢ = xᵢ - xᵢ₋₁
    c₀ = fᵢ₋₁
    c₁ = dfᵢ₋₁
    c₂ = @. (3 * (fᵢ - fᵢ₋₁) / hᵢ - (dfᵢ + 2 * dfᵢ₋₁)) / hᵢ
    c₃ = @. ((dfᵢ + dfᵢ₋₁) - 2 * (fᵢ - fᵢ₋₁) / hᵢ) / hᵢ^2
    s = @. c₀ + c₁ * (x - xᵢ₋₁) + c₂ * (x - xᵢ₋₁)^2 + c₃ * (x - xᵢ₋₁)^3
    return s
end

"""
    (solution::RungeKuttaSolution)(tₚ::Real)

interpolates `solution.u` at `tₚ` using linear splines.
"""
function (solution::RungeKuttaSolution)(tₚ::Real)
    @↓ u, t = solution
    n = findindex(tₚ, t, u)
    s = linearspline(tₚ, t[n-1], t[n], u[n-1], u[n])
    return s
end

"""
    (solution::RungeKuttaSolution)(tₚ::Real, f::Function)

interpolates `solution.u` at `tₚ` using Hermite's cubic splines.
"""
function (solution::RungeKuttaSolution)(tₚ::Real, f::Function)
    @↓ u, t = solution
    n = findindex(tₚ, t, u)
    duₙ₋₁, duₙ = f(u[n-1], t[n-1]), f(u[n], t[n])
    s = hermitecubicspline(tₚ, t[n-1], t[n], u[n-1], u[n], duₙ₋₁, duₙ)
    return s
end

#####
##### Functions
#####

"""
    length(solution::RungeKuttaSolution)

returns the number of time steps of `solution`.
"""
function Base.length(solution::RungeKuttaSolution)
    @↓ t = solution
    return length(t)
end

"""
    size(solution::RungeKuttaSolution)

returns a tuple containing the number of variables and of time steps of `solution`.
"""
function Base.size(solution::RungeKuttaSolution)
    @↓ u, t = solution
    return (length(u[1]), length(t))
end

"""
    extract(solution::RungeKuttaSolution, i::Integer) :: RungeKuttaSolution

returns the `i`th variable of `solution` (`i = 0` returns `t`).
"""
function extract(solution::RungeKuttaSolution, i::Integer)
    N = length(solution)
    @↓ u, t = solution
    if i == 0
        return t
    else
        return [u[n][i] for n = 1:N]
    end
end

"""
    extract(solution::RungeKuttaSolution, v::AbstractVector) :: RungeKuttaSolution

returns the variables of `solution` indicated in `v`.
"""
extract(solution::RungeKuttaSolution, v::AbstractVector) = tuple([extract(solution, i) for i in v]...)

"""
    extract(solution::RungeKuttaSolution) :: RungeKuttaSolution

returns all variables of `solution`, `t` included.
"""
function extract(solution::RungeKuttaSolution)
    d = size(solution)[1]
    return extract(solution, 0:d)
end

"""
    getindex(solution::RungeKuttaSolution, i::Integer) :: RungeKuttaSolution

returns a [`RungeKuttaSolution`](@ref) containing the fields of `solution` indexed at `i`.
"""
function Base.getindex(solution::RungeKuttaSolution, i::Integer)
    @↓ u, t = solution
    RungeKuttaSolution(u[i], t[i])
end

"""
    getindex(solution::RungeKuttaSolution, v::AbstractVector) :: RungeKuttaSolution

returns a [`RungeKuttaSolution`](@ref) containing the fields of `solution` indexed at `v`.
"""
function Base.getindex(solution::RungeKuttaSolution, v::AbstractVector)
    @↓ u, t = solution
    return RungeKuttaSolution(u[v], t[v])
end

"""
    setindex!(solution::RungeKuttaSolution, tuple::Tuple, i::Integer)

stores the values of `tuple` into the fields of `solution` indexed at `i`.
"""
function Base.setindex!(solution::RungeKuttaSolution, tuple::Tuple, i::Integer)
    @↓ u, t = solution
    ʊ, τ = tuple
    u[i] = ʊ
    t[i] = τ
    return solution
end

"""
    setindex!(solution::RungeKuttaSolution, point::RungeKuttaSolution, i::Integer)

stores the fields of `point` into the fields of `solution` indexed at `i`.
"""
function Base.setindex!(solution::RungeKuttaSolution, point::RungeKuttaSolution, i::Integer)
    @↓ u, t = solution
    @↓ ʊ, τ = point
    u[i] = ʊ
    t[i] = τ
    return solution
end

"""
    setindex!(solution::RungeKuttaSolution, smallersolution::RungeKuttaSolution, v::AbstractVector)

stores the fields of `smallersolution` into the fields of `solution` indexed at `v`.
"""
function Base.setindex!(solution::RungeKuttaSolution, smallersolution::RungeKuttaSolution, v::AbstractVector)
    @↓ u, t = solution
    @↓ ʊ, τ = smallersolution
    @. u[v] = ʊ
    @. t[v] = τ
    return solution
end

"""
    firstindex(solution::RungeKuttaSolution)

returns the first index of `solution`.
"""
function Base.firstindex(solution::RungeKuttaSolution)
    @↓ t = solution
    return firstindex(t)
end

"""
    lastindex(solution::RungeKuttaSolution)

returns the last index of `solution`.
"""
function Base.lastindex(solution::RungeKuttaSolution)
    @↓ t = solution
    return lastindex(t)
end
