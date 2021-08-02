"""
    RungeKuttaSolution{u_T, t_T, k_T} <: InitialValueSolution

returns a constructor for the numerical solution of an `InitialValueProblem`.

---

    RungeKuttaSolution(u, t[, k])

returns a `RungeKuttaSolution` with:
- `u` : numerical solution.
- `t` : time grid.
- `k` : stages.

---

    RungeKuttaSolution(problem::InitialValueProblem, solver::RungeKuttaSolver; save_stages::Bool = false)

returns an initialised `RungeKuttaSolution` given an `InitialValueProblem` and a `RungeKuttaSolver`. `save_stages` flags when to save all stages into `k`.
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

Base.summary(io::IO, solution::RungeKuttaSolution) = print(io, "RungeKuttaSolution")

function Base.show(io::IO, solution::RungeKuttaSolution)
    print(io, "RungeKuttaSolution:\n")
    pad = get(io, :pad, "")
    newline = get(io, :newline, "\n")
    names = propertynames(solution)
    N = length(names)
    for (n, name) in enumerate(names)
        field = getproperty(solution, name)
        print(io, pad, "   ‣ " * string(name) * " := ")
        show(IOContext(io, :pad => "   ", :newline => ""), field)
        n == N ? print(io, newline) : print(io, "\n")
    end
end

Base.length(solution::RungeKuttaSolution) = length(solution.t)
Base.getindex(solution::RungeKuttaSolution, k::Int) = RungeKuttaSolution(solution.u[k], solution.t[k])
Base.getindex(solution::RungeKuttaSolution, K::Vararg{Int, N}) where N = RungeKuttaSolution(solution.u[K], solution.t[K])
Base.getindex(solution::RungeKuttaSolution, K...) = RungeKuttaSolution(solution.u[K...], solution.t[K...])
function Base.setindex!(solution::RungeKuttaSolution, value::Tuple, k::Int)
    solution.u[k] = value[1]
    solution.t[k] = value[2]
    return solution
end
Base.lastindex(solution::RungeKuttaSolution) = lastindex(solution.t)

function extract(solution::RungeKuttaSolution, i)
    @↓ u, t = solution
    L = length(u[1])
    N = length(t)
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
