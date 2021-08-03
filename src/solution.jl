"""
    RungeKuttaSolution <: InitialValueSolution

A composite type for the solution of an [`InitialValueProblem`](@ref) obtained with a [`RungeKuttaSolver`](@ref).

# Constructors
```julia
RungeKuttaSolution(u, t[, k])
RungeKuttaSolution(problem, solver; save_stages = false)
```

# Arguments
- `u` : numerical solution.
- `t` : time grid.
- `k` : stages.
- `problem :: InitialValueProblem`.
- `solver :: RungeKuttaSolver`.
- `save_stages :: Bool` : flags when to save all stages into `k`.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct RungeKuttaSolution{u_T, t_T, k_T} <: InitialValueSolution
    u::u_T
    t::t_T
    k::k_T
end

RungeKuttaSolution(u, t) = RungeKuttaSolution(u, t, nothing)

function RungeKuttaSolution(problem::InitialValueProblem, solver::RungeKuttaSolver; save_stages::Bool = false)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ h = solver.stepsize
    @↓ s = solver.tableau
    N = round(Int, (tN - t0) / h) + 1
    u0_T = eltype(u0)
    L = length(u0)
    u = Vector{u0_T}(undef, N, L); u[1] = u0
    t = Vector{typeof(t0)}(undef, N); t[1] = t0
    k = save_stages ? Vector{u0_T}(undef, N, s, L) : nothing
    return RungeKuttaSolution(u, t, k)
end

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

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
    getindex(solution::RungeKuttaSolution, i::Integer)

returns a [`RungeKuttaSolution`](@ref) containing the fields of `solution` indexed at `i`.
"""
function Base.getindex(solution::RungeKuttaSolution, i::Integer)
    @↓ u, t, k = solution
    if k isa Nothing
        return RungeKuttaSolution(u[i], t[i])
    else
        return RungeKuttaSolution(u[i], t[i], k[i])
    end
end

"""
    setindex!(solution::RungeKuttaSolution, tuple::Tuple, i::Integer)

stores the values of `tuple` into the fields of `solution` indexed at `i`.
"""
function Base.setindex!(solution::RungeKuttaSolution, tuple::Tuple, i::Integer)
    @↓ u, t, k = solution
    u[i] = tuple[1]
    t[i] = tuple[2]
    if length(tuple) > 2 && !(k isa Nothing)
        k[i] = tuple[3]
    end
    return solution
end

Base.lastindex(solution::RungeKuttaSolution) = lastindex(solution.t)

"""
    show(io::IO, solution::RungeKuttaSolution)

prints a full description of `solution` and its contents to a stream `io`.
"""
Base.show(io::IO, solution::RungeKuttaSolution) = _show(io, solution)

"""
    summary(io::IO, solution::RungeKuttaSolution)

prints a brief description of `solution` to a stream `io`.
"""
Base.summary(io::IO, solution::RungeKuttaSolution) = _summary(io, solution)

# ---------------------------------------------------------------------------- #
#                                    Methods                                   #
# ---------------------------------------------------------------------------- #

function extract(solution::RungeKuttaSolution, i)
    @↓ u, t = solution
    (L, N) = size(solution)
    return RungeKuttaSolution([u[n][i] for n = 1:N], t)
end

# Make solution callable
function find_closest(target, vector)
    N = length(vector)
    if target ≤ vector[1]
        return (1, vector[1])
    elseif target ≥ vector[end]
        return (N, vector[end])
    end
    i = 0; j = N; n = 0
    while i < j
        n = (i + j) ÷ 2
        if target == vector[n]
            return vector[n]
        elseif target < vector[n]
            if (n > 1) && (target > vector[n - 1])
                if target - vector[n - 1] ≥ vector[n] - target
                    return (n, vector[n])
                else
                    return (n - 1, vector[n - 1])
                end
            end
            j = n
        else
            if (n < N) && (target < vector[n + 1])
                if target - vector[n] ≥ vector[n + 1] - target
                    return (n + 1, vector[n + 1])
                else
                    return (n, vector[n])
                end
            end
            i = n + 1
        end
    end
    return n, vector[n]
end

function (solution::RungeKuttaSolution)(t::Real)
    (n, tₙ) = find_closest(t, solution.t)
    uₙ = solution.u[n]
    return tₙ, uₙ
end
