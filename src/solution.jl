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
- [`dimension`](@ref) : number of variables.
"""
struct RungeKuttaSolution{u_T<:(AbstractVector{𝕍} where 𝕍<:AbstractVector{ℂ} where ℂ<:Number), t_T<:(AbstractVector{ℝ} where ℝ<:Real)} <: AbstractRungeKuttaSolution
    u::u_T
    t::t_T
end

function RungeKuttaSolution(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ h = solver.stepsize
    N = ceil(Int, (tN - t0) / h) + 1 # e.g. tspan = (0, 1), h = 0.3 ⇒ t = [0.0, 0.3, 0.6, 0.9, 1.2]
    u = [similar(u0) for i = 1:N]; u[1] = u0
    t = Vector{typeof(t0)}(undef, N); t[1] = t0
    return RungeKuttaSolution(u, t)
end

#----------------------------------- METHODS -----------------------------------

"""
    (solution::RungeKuttaSolution)(tₚ::Real)

interpolates `solution` using linear splines, approximating its value at `tₚ`.
"""
function (solution::RungeKuttaSolution)(tₚ::Real)
    @↓ u, t = solution
    N = length(t)
    if tₚ < t[1]
        return u[1]
    elseif tₚ ≥ t[N]
        return u[N]
    else
        for n in 2:N
            if t[n-1] ≤ tₚ < t[n]
                uₚ = linearspline(tₚ, t[n-1], t[n], u[n-1], u[n])
                return uₚ
            end
        end
    end
end

"""
    (solution::RungeKuttaSolution)(tₚ::Real, f::Function)

interpolates `solution` using Hermite's cubic splines, approximating its value at `tₚ`.
"""
function (solution::RungeKuttaSolution)(tₚ::Real, f::Function)
    @↓ u, t = solution
    N = length(t)
    if tₚ < t[1]
        return u[1]
    elseif tₚ ≥ t[N]
        return u[N]
    else
        for n in 2:N
            if t[n-1] ≤ tₚ < t[n]
                duₙ₋₁, duₙ = f(u[n-1], t[n-1]), f(u[n], t[n])
                uₚ = hermitecubicspline(tₚ, t[n-1], t[n], u[n-1], u[n], duₙ₋₁, duₙ)
                return uₚ
            end
        end
    end
end

#---------------------------------- FUNCTIONS ----------------------------------

"""
    length(solution::RungeKuttaSolution)

returns the number of time steps of `solution`.
"""
function Base.length(solution::RungeKuttaSolution)
    @↓ t = solution
    return length(t)
end

"""
    dimension(solution::RungeKuttaSolution)

returns the number of variables of `solution`.
"""
function dimension(solution::RungeKuttaSolution)
    @↓ u = solution
    return length(u[1])
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
    d = dimension(solution)
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
