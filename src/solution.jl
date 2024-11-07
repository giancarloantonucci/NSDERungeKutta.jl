"""
    RungeKuttaSolution <: AbstractRungeKuttaSolution

A composite type for an [`AbstractRungeKuttaSolution`](@ref) obtained using an [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
RungeKuttaSolution(u, t)
RungeKuttaSolution(problem, solver)
```

# Arguments
- `u :: AbstractVector{<:AbstractVector{<:Number}}` : numerical solution
- `t :: AbstractVector{<:Real}` : time grid

# Functions
- [`extract`](@ref) : extract all values for a specific variable
- [`firstindex`](@ref) : get the first index
- [`getindex`](@ref) : get specified value(s) and time
- [`lastindex`](@ref) : get the last index
- [`length`](@ref) : get the number of time steps
- [`setindex!`](@ref) : set value(s) and time
- [`numtimesteps`](@ref) : get the number of time steps
- [`numvariables`](@ref) : get the number of variables
"""
struct RungeKuttaSolution{u_T<:AbstractVector{<:AbstractVector{<:Number}}, t_T<:AbstractVector{<:Real}} <: AbstractRungeKuttaSolution
    u::u_T
    t::t_T
end

function RungeKuttaSolution(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ h = solver.stepsize
    N = ceil(Int, (tN - t0) / h) + 1 # e.g. tspan = (0, 1), h = 0.3 -> t = [0.0, 0.3, 0.6, 0.9, 1.2]
    u = [similar(u0) for _ = 1:N]
    u[1] = u0
    t = Vector{typeof(t0)}(undef, N)
    t[1] = t0
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

uses Hermite's cubic splines to interpolate `solution` and approximate its value at `tₚ`. Note that it needs the derivative function `f(u, t)`, e.g. from an [`NSDEBase.AbstractRightHandSide`](@extref) subtype.
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

returns the number of time steps in `solution`.
"""
Base.length(solution::RungeKuttaSolution) = length(solution.t)

"""
    numtimesteps(solution::RungeKuttaSolution)

returns the number of time steps in `solution`.
"""
numtimesteps(solution::RungeKuttaSolution) = length(solution.t)

"""
    numvariables(solution::RungeKuttaSolution)

returns the number of variables in `solution`.
"""
numvariables(solution::RungeKuttaSolution) = length(solution.u[1])

"""
    size(solution::RungeKuttaSolution)

returns a tuple containing the number of variables and time steps in `solution`.
"""
Base.size(solution::RungeKuttaSolution) = (numvariables(solution), numtimesteps(solution))

"""
    extract(solution::RungeKuttaSolution, i::Integer) :: RungeKuttaSolution

returns the `i`-th variable of `solution`. `i = 0` returns `t`.
"""
extract(solution::RungeKuttaSolution, i::Integer) = i == 0 ? solution.t : [solution.u[n][i] for n = 1:length(solution)]

"""
    extract(solution::RungeKuttaSolution, v::AbstractVector) :: RungeKuttaSolution

returns the variables of `solution` indicated by the indices `v`.
"""
extract(solution::RungeKuttaSolution, v::AbstractVector) = tuple([extract(solution, i) for i in v]...)

"""
    extract(solution::RungeKuttaSolution) :: RungeKuttaSolution

returns all variables of `solution`, including `t`.
"""
extract(solution::RungeKuttaSolution) = extract(solution, 0:numvariables(solution))

"""
    getindex(solution::RungeKuttaSolution, i::Integer) :: RungeKuttaSolution

returns new a [`RungeKuttaSolution`](@ref) containing the fields of `solution` at index `i`.
"""
Base.getindex(solution::RungeKuttaSolution, i::Integer) = RungeKuttaSolution(solution.u[i], solution.t[i])

"""
    getindex(solution::RungeKuttaSolution, v::AbstractVector) :: RungeKuttaSolution

returns a new [`RungeKuttaSolution`](@ref) containing the fields of `solution` at the indices `v`.
"""
Base.getindex(solution::RungeKuttaSolution, v::AbstractVector) = RungeKuttaSolution(solution.u[v], solution.t[v])

"""
    setindex!(solution::RungeKuttaSolution, values::Tuple, i::Integer)

stores the values from `values` into the fields of `solution` at the specified index `i`.
"""
function Base.setindex!(solution::RungeKuttaSolution, values::Tuple, i::Integer)
    @↓ u, t = solution
    u_new, t_new = values
    u[i] = u_new
    t[i] = t_new
    return solution
end

"""
    setindex!(solution::RungeKuttaSolution, values::RungeKuttaSolution, i::Integer)

stores the fields of `values` into the fields of `solution` at the specified index `i`.
"""
function Base.setindex!(solution::RungeKuttaSolution, values::RungeKuttaSolution, i::Integer)
    @↓ u, t = solution
    @↓ u_new, t_new = values
    u[i] = u_new
    t[i] = t_new
    return solution
end

"""
    setindex!(solution::RungeKuttaSolution, values::RungeKuttaSolution, v::AbstractVector)

stores the fields of `values` into the fields of `solution` at the specified indices `v`.
"""
function Base.setindex!(solution::RungeKuttaSolution, values::RungeKuttaSolution, v::AbstractVector)
    @↓ u, t = solution
    @↓ u_new, t_new = values
    @. u[v] = u_new
    @. t[v] = t_new
    return solution
end

"""
    firstindex(solution::RungeKuttaSolution)

returns the first index of `solution`.
"""
Base.firstindex(solution::RungeKuttaSolution) = firstindex(solution.t)

"""
    lastindex(solution::RungeKuttaSolution)

returns the last index of `solution`.
"""
Base.lastindex(solution::RungeKuttaSolution) = lastindex(solution.t)
