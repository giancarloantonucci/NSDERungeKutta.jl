"""
    RungeKuttaSolution(u, t[, k]) -> RungeKuttaSolution

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
    u = Vector{eltype(u0)}(undef, N, length(u0))
    u[1] = u0
    t = Vector{typeof(t0)}(undef, N)
    t[1] = t0
    k = save_stages ? Vector{uT}(undef, N, s, L) : nothing
    return RungeKuttaSolution(u, t, k)
end

Base.length(solution::RungeKuttaSolution) = length(solution.t)

function findclosest(t, target)
    N = length(t)
    if target ≤ t[1]
        return (1, t[1])
    elseif target ≥ t[end]
        return (N, t[end])
    end
    i = 0; j = N; n = 0
    while i < j
        n = (i + j) ÷ 2
        if target == t[n]
            return t[n]
        elseif target < t[n]
            if (n > 1) && (target > t[n - 1])
                if target - t[n - 1] ≥ t[n] - target
                    return (n, t[n])
                else
                    return (n - 1, t[n - 1])
                end
            end
            j = n
        else
            if (n < N) && (target < t[n + 1])
                if target - t[n] ≥ t[n + 1] - target
                    return (n + 1, t[n + 1])
                else
                    return (n, t[n])
                end
            end
            i = n + 1
        end
    end
    return (n, t[n])
end

function (solution::RungeKuttaSolution)(t::Real)
    (n, tₙ) = findclosest(solution.t, t)
    uₙ = solution.u[n]
    return tₙ, uₙ
end
