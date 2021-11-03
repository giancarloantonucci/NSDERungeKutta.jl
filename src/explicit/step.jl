"""
    step!(cache::ExplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver)

computes a step of the [`AbstractRungeKuttaSolution`](@ref) of an [`AbstractInitialValueProblem`](@ref) using an [`ExplicitRungeKuttaSolver`](@ref).
"""
function step!(cache::ExplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, problem::AbstractInitialValueProblem, solver::ExplicitRungeKuttaSolver)
    @↓ n, v = cache
    @↓ u, t, savestages = solution
    k = savestages ? solution.k[n] : cache.k
    @↓ f! = problem.rhs
    @↓ s, A, c, b = solver.tableau
    @↓ h = solver.stepsize
    v = u[n+1] # avoid allocs
    # Compute stages
    for i = 1:s
        zero!(v)
        for j = 1:(i-1)
            @. v += A[i,j] * k[j]
        end
        @. v = u[n] + h * v
        # TO-DO: @← k[i] = f(v, t[n] + h * c[i])
        f!(k[i], v, t[n] + h * c[i])
    end
    # Compute step
    zero!(v)
    for i = 1:s
        @. v += b[i] * k[i]
    end
    @. v = u[n] + h * v
    t[n+1] = t[n] + h
    return u[n+1], t[n+1]
end
