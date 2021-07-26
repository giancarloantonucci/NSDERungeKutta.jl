"""
    RungeKuttaSolution(u, t[, k]) <: InitialValueSolution

returns a constructor for the numerical solution of an `InitialValueProblem`.

# Arguments
- `u` : numerical solution.
- `t` : time grid.
- `k` : stages.
"""
struct RungeKuttaSolution{u_T, t_T, k_T} <: InitialValueSolution
    u::u_T
    t::t_T
    k::k_T
end

RungeKuttaSolution(u, t) = RungeKuttaSolution(u, t, nothing)

function RungeKuttaSolution(problem::InitialValueProblem, solver::RungeKuttaSolver, save_stages)
    @↓ u0, (t0, tN) ← tspan = problem
    @↓ h = solver.stepsize
    @↓ s = solver.tableau
    N = floor(Int, (tN - t0) / h) + 1
    u0_T = eltype(u0)
    u = Vector{u0_T}(undef, N, length(u0))
    u[1] = u0
    t = Vector{typeof(t0)}(undef, N)
    t[1] = t0
    k = if save_stages
        Vector{u0_T}(undef, N, s, L)
    else
        nothing
    end
    return RungeKuttaSolution(u, t, k)
end

# Utility functions
Base.length(solution::RungeKuttaSolution) = length(solution.t)

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
