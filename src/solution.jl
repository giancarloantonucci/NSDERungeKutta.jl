# NSDERungeKutta/src/solution.jl

"""
    RungeKuttaSolution <: AbstractRungeKuttaSolution

A composite type for an [`AbstractRungeKuttaSolution`](@ref) obtained using an [`AbstractRungeKuttaSolver`](@ref).

# Constructors
```julia
RungeKuttaSolution(u, t, k)
RungeKuttaSolution(problem, solver; dense=false)
```

# Arguments
- `u :: AbstractVector{<:AbstractVector{<:Number}}` : numerical solution
- `t :: AbstractVector{<:Real}` : time grid
- `k :: AbstractVector{<:AbstractVector{<:AbstractVector{<:Number}}}` : stages history (for dense output)

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
struct RungeKuttaSolution{
            u_T <: AbstractVector{<:AbstractVector{<:Number}},
            t_T <: AbstractVector{<:Real},
            k_T <: Union{AbstractVector{<:AbstractVector{<:AbstractVector{<:Number}}}, Nothing}
        } <: AbstractRungeKuttaSolution
    u :: u_T
    t :: t_T
    k :: k_T
end

function RungeKuttaSolution(problem::AbstractInitialValueProblem, solver::AbstractRungeKuttaSolver; dense::Bool=false)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ h = solver.stepsize
    N = ceil(Int, (tN - t0) / h) + 1 # e.g. tspan = (0, 1), h = 0.3 -> t = [0.0, 0.3, 0.6, 0.9, 1.2]
    u = [similar(u0) for _ = 1:N]
    copyto!(u[1], u0)
    t = Vector{typeof(t0)}(undef, N)
    t[1] = t0

    if dense
        # Outer: time, Middle: stages, Inner: state
        k = Vector{Vector{typeof(u0)}}()
        sizehint!(k, N)
        return RungeKuttaSolution(u, t, k)
    else
        return RungeKuttaSolution(u, t, nothing)
    end
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

"""
    (solution::RungeKuttaSolution)(tₚ::Real, tableau::AbstractButcherTableau)

Evaluates the dense output solution at `tₚ` using the stored stages and tableau coefficients.
"""
function (solution::RungeKuttaSolution)(tₚ::Real, tableau::AbstractButcherTableau)
    @↓ u, t, k = solution
    
    if k isa Nothing || tableau.b_dense isa Nothing
        # Fallback to linear spline if data is missing
        return solution(tₚ) 
    end

    N = length(t)
    # Find index n such that t[n] <= tₚ <= t[n+1]
    n = searchsortedlast(t, tₚ)
    
    if n == N
        return u[N]
    end

    if n == 0
        return u[1]
    end

    # Calculation for dense output
    dt = t[n+1] - t[n]
    θ = (tₚ - t[n]) / dt
    
    # b(θ) calculation
    # b_dense is a matrix where column j contains coeffs for stage j
    # b_j(θ) = b_dense[1,j]*θ + b_dense[2,j]*θ^2 + ...
    
    uₚ = copy(u[n]) # Start with u_n
    
    # Perform the summation: u(θ) = u_n + h * Σ b_j(θ) * k_j
    # Note: k[n] contains the stages for step n
    
    stages = k[n]
    s = length(stages)
    
    for j = 1:s
        # Evaluate polynomial for weight j
        # coeffs = tableau.b_dense[:, j]
        bj_θ = 0.0
        θ_pow = θ
        for p in axes(tableau.b_dense, 1)
             bj_θ += tableau.b_dense[p, j] * θ_pow
             θ_pow *= θ
        end
        
        @. uₚ += dt * bj_θ * stages[j]
    end
    
    return uₚ
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
function Base.getindex(solution::RungeKuttaSolution, i::Integer)
    @↓ u, t, k = solution
    new_u = [u[i]]
    new_t = [t[i]]
    new_k = k isa Nothing ? nothing : [k[i]]
    return RungeKuttaSolution(new_u, new_t, new_k)
end

"""
    getindex(solution::RungeKuttaSolution, v::AbstractVector) :: RungeKuttaSolution

returns a new [`RungeKuttaSolution`](@ref) containing the fields of `solution` at the indices `v`.
"""
function Base.getindex(solution::RungeKuttaSolution, v::AbstractVector)
    @↓ u, t, k = solution
    new_k = k isa Nothing ? nothing : k[v]
    return RungeKuttaSolution(solution.u[v], solution.t[v], new_k)
end

"""
    setindex!(solution::RungeKuttaSolution, values::Tuple, i::Integer)

stores the values from `values` into the fields of `solution` at the specified index `i`.
If `solution` is dense (has `k`), `values` must be a 3-tuple `(u, t, k)`.
Otherwise, `values` must be a 2-tuple `(u, t)`.
"""
function Base.setindex!(solution::RungeKuttaSolution, values::Tuple, i::Integer)
    @↓ u, t, k = solution
    if k isa Nothing
        u_new, t_new = values
        u[i] = u_new
        t[i] = t_new
    else
        u_new, t_new, k_new = values
        u[i] = u_new
        t[i] = t_new
        k[i] = k_new
    end
    return solution
end

"""
    setindex!(solution::RungeKuttaSolution, values::RungeKuttaSolution, i::Integer)

stores the fields of `values` into the fields of `solution` at the specified index `i`.
Assumes `values` contains data for a single time step (e.g. from `getindex`).
"""
function Base.setindex!(solution::RungeKuttaSolution, values::RungeKuttaSolution, i::Integer)
    @↓ u, t, k = solution
    @↓ u_new, t_new, k_new = values
    u[i] = u_new
    t[i] = t_new
    if !(k isa Nothing)
        if k_new isa Nothing
            error("Cannot assign non-dense solution data to a dense solution at index $i.")
        else
            k[i] = k_new
        end
    end
    return solution
end

"""
    setindex!(solution::RungeKuttaSolution, values::RungeKuttaSolution, v::AbstractVector)

stores the fields of `values` into the fields of `solution` at the specified indices `v`.
"""
function Base.setindex!(solution::RungeKuttaSolution, values::RungeKuttaSolution, v::AbstractVector)
    @↓ u, t, k = solution
    @↓ u_new, t_new, k_new = values
    @. u[v] = u_new
    @. t[v] = t_new
    if !(k isa Nothing)
        if k_new isa Nothing
            error("Cannot assign non-dense solution data to a dense solution.")
        else
            @. k[v] = k_new
        end
    end
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
