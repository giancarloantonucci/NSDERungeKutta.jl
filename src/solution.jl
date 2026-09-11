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
    end
    # Binary search, like the dense-output path: this call sits inside the
    # Parareal hot loop (chunk boundary values, every sweep), where a linear
    # scan over a 10³-step chunk multiplies through steps × s × N × K × M.
    n = searchsortedlast(t, tₚ) # t[n] ≤ tₚ < t[n+1]
    return linearspline(tₚ, t[n], t[n+1], u[n], u[n+1])
end

"""
    (solution::RungeKuttaSolution)(tₚ::Real, f::Function)

uses Hermite's cubic splines to interpolate `solution` and approximate its value at `tₚ`. Note that it needs the derivative function `f(u, t)`, e.g. from an `NSDEBase.AbstractRightHandSide` subtype.
"""
function (solution::RungeKuttaSolution)(tₚ::Real, f::Function)
    @↓ u, t = solution
    N = length(t)
    if tₚ < t[1]
        return u[1]
    elseif tₚ ≥ t[N]
        return u[N]
    end
    n = searchsortedlast(t, tₚ) # t[n] ≤ tₚ < t[n+1]
    duₙ, duₙ₊₁ = f(u[n], t[n]), f(u[n+1], t[n+1])
    return hermitecubicspline(tₚ, t[n], t[n+1], u[n], u[n+1], duₙ, duₙ₊₁)
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

# Stage history lives on INTERVALS, not nodes: `k[n]` holds the stages of the
# step from `t[n]` to `t[n+1]`, so a dense solution with N nodes has N − 1
# entries in `k` (see `solve!`). Every indexer below honours that: a single node
# owns no interval, a contiguous range of m nodes owns m − 1 intervals, and a
# non-contiguous selection owns none that could be reused (the polynomials
# between non-adjacent nodes are not the ones stored), so it is refused rather
# than fitted with the wrong stages.

"""
    getindex(solution::RungeKuttaSolution, i::Integer) :: RungeKuttaSolution

returns a new [`RungeKuttaSolution`](@ref) containing the fields of `solution`
at node `i`. A single node carries no interval, so the slice has no dense stage
history even when `solution` does.
"""
function Base.getindex(solution::RungeKuttaSolution, i::Integer)
    @↓ u, t = solution
    return RungeKuttaSolution([u[i]], [t[i]], nothing)
end

"""
    getindex(solution::RungeKuttaSolution, v::AbstractUnitRange) :: RungeKuttaSolution

returns a new [`RungeKuttaSolution`](@ref) containing the fields of `solution`
at the contiguous nodes `v`, together with the stage history of the
`length(v) − 1` intervals between them (when `solution` is dense).
"""
function Base.getindex(solution::RungeKuttaSolution, v::AbstractUnitRange)
    @↓ u, t, k = solution
    new_k = (k isa Nothing || isempty(v)) ? nothing : k[first(v):last(v)-1]
    return RungeKuttaSolution(u[v], t[v], new_k)
end

"""
    getindex(solution::RungeKuttaSolution, v::AbstractVector) :: RungeKuttaSolution

returns a new [`RungeKuttaSolution`](@ref) containing the fields of `solution`
at the nodes `v`. For a non-contiguous selection of a dense solution the
stored stage polynomials do not describe the gaps between the chosen nodes,
so the request is refused with an `ArgumentError`; slice `solution.u` and
`solution.t` directly if only the nodes are wanted.
"""
function Base.getindex(solution::RungeKuttaSolution, v::AbstractVector)
    @↓ u, t, k = solution
    if !(k isa Nothing) && !(isempty(v) || v == first(v):last(v))
        throw(ArgumentError("cannot slice a dense `RungeKuttaSolution` at non-contiguous nodes: the stored stages belong to adjacent intervals only. Use a range, or index `solution.u`/`solution.t` directly."))
    end
    new_k = (k isa Nothing || isempty(v)) ? nothing : k[first(v):last(v)-1]
    return RungeKuttaSolution(u[v], t[v], new_k)
end

# Writing through the indexer into a DENSE solution: the stage history is a
# per-interval record that must agree with the nodes on both sides of each
# interval. Changing a node without the two adjoining intervals' stages
# leaves those intervals describing a curve through the OLD node, and
# `solution(t, tableau)` would then silently interpolate off the stored
# trajectory. There is no way to repair that from node data alone, so the
# indexer refuses partial writes into a dense solution; either write the whole
# solution (all nodes with all `N − 1` interval stages), or edit `u`, `t`
# and `k` directly and take responsibility for their consistency. Non-dense
# solutions carry no such record and accept any node write.

_dense_partial_write() = throw(ArgumentError(
    "cannot write single nodes into a dense `RungeKuttaSolution`: the stored stages of the adjoining intervals would no longer match the nodes. Write the whole solution (all nodes and all interval stages), or edit `solution.u`, `solution.t` and `solution.k` directly."))

_contiguous(v::AbstractVector) = isempty(v) || v == first(v):last(v)

"""
    setindex!(solution::RungeKuttaSolution, values::Tuple, i::Integer)

stores `(u, t)` from `values` at node `i` of a non-dense `solution`. A dense
solution refuses the write (see the note on partial writes above).
"""
function Base.setindex!(solution::RungeKuttaSolution, values::Tuple, i::Integer)
    @↓ u, t, k = solution
    k isa Nothing || _dense_partial_write()
    length(values) == 2 || throw(ArgumentError("expected a `(u, t)` pair, got a $(length(values))-tuple."))
    u_new, t_new = values
    u[i] = u_new
    t[i] = t_new
    return solution
end

"""
    setindex!(solution::RungeKuttaSolution, values::RungeKuttaSolution, i::Integer)

stores the single node held by `values` (e.g. from `getindex`) at node `i` of
a non-dense `solution`. A dense solution refuses the write (see the note on
partial writes above).
"""
function Base.setindex!(solution::RungeKuttaSolution, values::RungeKuttaSolution, i::Integer)
    @↓ u, t, k = solution
    k isa Nothing || _dense_partial_write()
    length(values) == 1 || throw(DimensionMismatch("expected a single-node slice, got $(length(values)) nodes."))
    u[i] = values.u[1]
    t[i] = values.t[1]
    return solution
end

"""
    setindex!(solution::RungeKuttaSolution, values::RungeKuttaSolution, v::AbstractVector)

stores the nodes of `values` at the nodes `v` of `solution`. On a non-dense
target any `v` is accepted. On a dense target the write must cover the WHOLE
solution — `v` equal to `1:length(solution)` and `values` dense with the
matching `length(v) − 1` interval stages — so that nodes and stages are
replaced together; anything less is refused (see the note above). All checks
run before anything is modified.
"""
function Base.setindex!(solution::RungeKuttaSolution, values::RungeKuttaSolution, v::AbstractVector)
    @↓ u, t, k = solution
    @↓ u_new ← u, t_new ← t, k_new ← k = values
    # Validate everything first: a refused write must leave `solution` intact.
    length(u_new) == length(v) && length(t_new) == length(v) ||
        throw(DimensionMismatch("expected $(length(v)) nodes, got $(length(u_new))."))
    all(i -> checkbounds(Bool, t, i), v) || throw(BoundsError(solution, v))
    if !(k isa Nothing)
        (_contiguous(v) && !isempty(v) && first(v) == firstindex(t) && last(v) == lastindex(t)) ||
            _dense_partial_write()
        k_new isa Nothing && throw(ArgumentError("a dense solution can only be overwritten by a dense one carrying its interval stages."))
        length(k_new) == length(v) - 1 ||
            throw(DimensionMismatch("expected $(length(v) - 1) stage entries for $(length(v)) nodes, got $(length(k_new))."))
    end
    for (j, i) in enumerate(v)
        u[i] = u_new[j]
        t[i] = t_new[j]
    end
    if !(k isa Nothing)
        intervals = first(v):last(v)-1 # computed OUTSIDE any broadcast: `@.` would dot `first`/`last`
        for (j, i) in enumerate(intervals)
            k[i] = k_new[j]
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
