"""
    RungeKuttaSolution <: AbstractRungeKuttaSolution

A composite type for an [`AbstractRungeKuttaSolution`](@ref).

# Constructors
```julia
RungeKuttaSolution(u, t[, k])
RungeKuttaSolution(problem, solver; savestages = false)
```

# Arguments
- `u          :: AbstractVector`              : numerical solution.
- `t          :: AbstractVector`              : time grid.
- `k          :: AbstractVector`              : vector of stages.
- `problem    :: AbstractInitialValueProblem` : initial value problem.
- `solver     :: AbstractRungeKuttaSolver`    : Runge-Kutta solver.
- `savestages :: Bool`                        : flags when to save all stages into `k`.

# Functions
- [`getindex` ](@ref) : get value(s) and time.
- [`lastindex`](@ref) : last index.
- [`length`   ](@ref) : number of time steps.
- [`setindex!`](@ref) : set value(s) and time.
- [`show`     ](@ref) : shows name and contents.
- [`size`     ](@ref) : number of variables and time steps.
- [`summary`  ](@ref) : shows name.
"""
mutable struct RungeKuttaSolution{u_T, t_T, k_T} <: AbstractRungeKuttaSolution
    u::u_T
    t::t_T
    k::k_T
end

RungeKuttaSolution(u, t) = RungeKuttaSolution(u, t, nothing)

function RungeKuttaSolution(
    problem::AbstractInitialValueProblem,
    solver::AbstractRungeKuttaSolver;
    savestages::Bool = false
)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ h = solver.stepsize
    N = round(Int, (tN - t0) / h) + 1
    # Initialise u
    u0_T = eltype(u0)
    L = length(u0)
    u = Vector{u0_T}(undef, N, L)
    u[1] = u0
    # Initialise t
    t0_T = typeof(t0)
    t = Vector{t0_T}(undef, N)
    t[1] = t0
    # Initialise k
    if savestages
        @↓ s = solver.tableau
        k = Vector{u0_T}(undef, N, s, L)
        return RungeKuttaSolution(u, t, k)
    else
        return RungeKuttaSolution(u, t)
    end
end

############################################################################################
#                                         FUNCTIONS                                        #
############################################################################################

"""
    length(solution::AbstractRungeKuttaSolution)

returns the number of time steps of `solution`.
"""
Base.length(solution::AbstractRungeKuttaSolution) = length(solution.t)

"""
    size(solution::AbstractRungeKuttaSolution)

returns a tuple containing the number of variables and of time steps of `solution`.
"""
Base.size(solution::AbstractRungeKuttaSolution) = (length(solution.u[1]), length(solution.t))

"""
    getindex(solution::AbstractRungeKuttaSolution, i::Integer)

returns a [`RungeKuttaSolution`](@ref) containing the fields of `solution` indexed at `i`.
"""
function Base.getindex(solution::AbstractRungeKuttaSolution, i::Integer)
    @↓ u, t, k = solution
    if k isa Nothing
        return RungeKuttaSolution(u[i], t[i])
    else
        return RungeKuttaSolution(u[i], t[i], k[i])
    end
end

"""
    setindex!(solution::AbstractRungeKuttaSolution, tuple::Tuple, i::Integer)

stores the values of `tuple` into the fields of `solution` indexed at `i`.
"""
function Base.setindex!(solution::AbstractRungeKuttaSolution, tuple::Tuple, i::Integer)
    @↓ u, t, k = solution
    u[i] = tuple[1]
    t[i] = tuple[2]
    if length(tuple) > 2 && !(k isa Nothing)
        k[i] = tuple[3]
    end
    return solution
end

"""
    lastindex(solution::AbstractRungeKuttaSolution)

returns the last index of `solution`.
"""
Base.lastindex(solution::AbstractRungeKuttaSolution) = lastindex(solution.t)

############################################################################################
#                                         PRINTING                                         #
############################################################################################

"""
    show(io::IO, solution::AbstractRungeKuttaSolution)

prints a full description of `solution` and its contents to a stream `io`.
"""
Base.show(io::IO, solution::AbstractRungeKuttaSolution) = NSDEBase._show(io, solution)

"""
    summary(io::IO, solution::AbstractRungeKuttaSolution)

prints a brief description of `solution` to a stream `io`.
"""
Base.summary(io::IO, solution::AbstractRungeKuttaSolution) = NSDEBase._summary(io, solution)

############################################################################################
#                                          METHODS                                         #
############################################################################################

function extract(solution::AbstractRungeKuttaSolution, i)
    (L, N) = size(solution)
    @↓ u, t = solution
    return RungeKuttaSolution([u[n][i] for n = 1:N], t)
end

# Make solution callable
function findclosest(target, vector)
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

function (solution::AbstractRungeKuttaSolution)(tₚ::Real)
    @↓ u, t = solution
    (n, tₙ) = findclosest(tₚ, t)
    return tₙ, u[n]
end

function (solution::AbstractRungeKuttaSolution)(tspan::Tuple{Real, Real})
    @↓ u, t = solution
    n₁, _ = findclosest(tspan[1], t)
    n₂, _ = findclosest(tspan[2], t)
    return RungeKuttaSolution(u[n₁:n₂], t[n₁:n₂])
end

(solution::AbstractRungeKuttaSolution)(t₁::Real, t₂::Real) = solution((t₁, t₂))
