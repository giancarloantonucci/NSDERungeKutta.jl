"A fixed schedule based on a nominal duration, independent of translated endpoints."
struct StepSchedule{T<:AbstractFloat}
    duration::T
    h::T
    full::Int
    tail::T
end
function StepSchedule(duration::Real, h::Real)
    dt, step = promote(float(duration), float(h))
    isfinite(dt) && dt >= 0 && isfinite(step) && step > 0 ||
        throw(ArgumentError("finite nonnegative duration and positive step required"))
    q = dt / step
    isfinite(q) && q < typemax(Int) ÷ 2 ||
        throw(ArgumentError("step count exceeds supported range"))
    nearest = round(Int, q)
    divides = nearest >= 1 && abs(q - nearest) <= 8eps(max(q, one(q)))
    n = divides ? nearest : floor(Int, q)
    tail = divides ? zero(dt) : dt - n * step
    # Fold a rounding-sized remainder into the preceding step. The last step
    # then spans the residual interval instead of merely relabelling its time.
    if n > 0 && 0 < tail <= 1e-6 * step
        n -= 1
        tail = dt - n * step
    end
    return StepSchedule(dt, step, n, tail)
end

"Prepared fixed-step chunk solver. Normally built through `chunk_solver`."
struct ScheduledChunk{S<:AbstractRungeKuttaSolver,P<:StepSchedule} <:
       NSDEBase.AbstractInitialValueSolver
    solver::S
    plan::P
end
function ScheduledChunk(solver::AbstractRungeKuttaSolver, duration::Real)
    duration > 0 || throw(ArgumentError("chunk duration must be positive"))
    solver.adaptive === nothing ||
        throw(ArgumentError("a fixed schedule requires a fixed-step method"))
    return ScheduledChunk(solver, StepSchedule(duration, solver.stepsize.h))
end
function NSDEBase.chunk_solver(solver::AbstractRungeKuttaSolver, duration::Real)
    return solver.adaptive === nothing ? ScheduledChunk(solver, duration) : solver
end
function NSDEBase.chunk_solver(solver::ScheduledChunk, duration::Real)
    return ScheduledChunk(solver.solver, duration)
end
# Mutable stage storage belongs to each solve, never to the shared schedule.
struct ScheduledCache{C,T} <: NSDEBase.AbstractInitialValueCache
    full::C
    tail::T
end
function tail_cache(p, solver, plan, full)
    return if solver isa ExponentialRungeKuttaSolver && plan.tail > 0
        if full.tail === nothing || full.tail.h != plan.tail
            full.tail = RungeKuttaCache(p, withstepsize(solver, plan.tail))
        end
        full.tail
    else
        full
    end
end
function NSDEBase.initialize_cache(p::NSDEBase.AbstractInitialValueProblem, c::ScheduledChunk)
    full = RungeKuttaCache(p, c.solver)
    tail = tail_cache(p, c.solver, c.plan, full)
    return ScheduledCache(full, tail)
end
function NSDEBase.initialize_solution(
    p::NSDEBase.AbstractInitialValueProblem, c::ScheduledChunk; dense::Bool=false
)
    return RungeKuttaSolution(p, c.solver; dense)
end

function ensure_capacity!(sol, n, u0)
    n <= length(sol.u) && return sol
    append!(sol.u, [similar(u0) for _ in (length(sol.u) + 1):n])
    resize!(sol.t, n)
    return sol
end
function reset_solution!(cache, sol, p)
    cache.n = 1
    cache.m = 1
    cache.e[] = 0
    ensure_capacity!(sol, 1, p.u0)
    copyto!(sol.u[1], p.u0)
    sol.t[1] = p.tspan[1]
    sol.k === nothing || empty!(sol.k)
    if sol.stepsizes !== nothing
        empty!(sol.stepsizes.accepted)
        empty!(sol.stepsizes.rejected)
        push!(sol.stepsizes.rejected, eltype(sol.stepsizes.accepted)[])
    end
    return sol
end
function finish_solution!(sol, n)
    resize!(sol.u, n)
    resize!(sol.t, n)
    return sol
end

function solve_fixed!(cache, sol, p, solver, plan, tail=tail_cache(p, solver, plan, cache))
    check_dense_support(solver, sol.k !== nothing)
    reset_solution!(cache, sol, p)
    M = plan.full + (plan.tail > 0)
    ensure_capacity!(sol, M + 1, p.u0)
    for j in 1:M
        shortened = j > plan.full
        hused = shortened ? plan.tail : plan.h
        work = shortened ? tail : cache
        work.n = j
        work.m = 1
        work.e[] = cache.e[]
        sol.t[j] + typeof(sol.t[j])(hused) > sol.t[j] ||
            throw(StepFailure("fixed step cannot advance time at this precision"))
        step!(work, sol, p, solver, hused)
        cache.e[] = work.e[]
        sol.k === nothing || push!(sol.k, copy.(work.k))
        if sol.stepsizes !== nothing
            push!(sol.stepsizes.accepted, hused)
            push!(sol.stepsizes.rejected, eltype(sol.stepsizes.accepted)[])
        end
    end
    cache.n = M + 1
    sol.t[M + 1] = p.tspan[2]
    return finish_solution!(sol, M + 1)
end
function NSDEBase.solve!(
    cache::ScheduledCache,
    sol::AbstractRungeKuttaSolution,
    p::NSDEBase.AbstractInitialValueProblem,
    c::ScheduledChunk,
)
    t0, tN = p.tspan
    scale = max(abs(float(t0)), abs(float(tN)), c.plan.duration)
    abs((tN - t0) - c.plan.duration) <= 32eps(scale) ||
        throw(ArgumentError("problem span differs from nominal chunk duration"))
    return solve_fixed!(cache.full, sol, p, c.solver, c.plan, cache.tail)
end
function NSDEBase.solve(
    p::NSDEBase.AbstractInitialValueProblem, c::ScheduledChunk; dense::Bool=false
)
    sol = NSDEBase.initialize_solution(p, c; dense)
    return NSDEBase.solve!(NSDEBase.initialize_cache(p, c), sol, p, c)
end
function (c::ScheduledChunk)(p::NSDEBase.AbstractInitialValueProblem; kwargs...)
    return NSDEBase.solve(p, c; kwargs...)
end
