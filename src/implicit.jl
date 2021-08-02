"""
    ImplicitRungeKuttaSolver{tableau_T, stepsize_T, newton_T, adaptive_T} <: RungeKuttaSolver

returns a constructor for an implicit `RungeKuttaSolver`.

---

    ImplicitRungeKuttaSolver(tableau, h, newton[, adaptive])
    IRK(args...; kwargs...) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` with:
- `tableau  :: ButcherTableau`     : Butcher tableau.
- `stepsize :: StepSize`           : step-size.
- `newton   :: NewtonParameters`   : simplified Newton's parameters
- `adaptive :: AdaptiveParameters` : embedded method's parameters.
"""
struct ImplicitRungeKuttaSolver{tableau_T, stepsize_T, newton_T, adaptive_T} <: RungeKuttaSolver
    tableau::tableau_T
    stepsize::stepsize_T
    newton::newton_T
    adaptive::adaptive_T
end

ImplicitRungeKuttaSolver(tableau, h::Real, newton, adaptive) = ImplicitRungeKuttaSolver(tableau, StepSize(h), newton, adaptive)
ImplicitRungeKuttaSolver(tableau, stepsize, newton) = ImplicitRungeKuttaSolver(tableau, stepsize, newton, nothing)
@doc (@doc ImplicitRungeKuttaSolver) IRK(args...; kwargs...) = ImplicitRungeKuttaSolver(args...; kwargs...)

Base.summary(io::IO, solver::ImplicitRungeKuttaSolver) = print(io, "ImplicitRungeKuttaSolver")

function Base.show(io::IO, solver::ImplicitRungeKuttaSolver)
    print(io, "ImplicitRungeKuttaSolver:")
    pad = get(io, :pad, "")
    names = propertynames(solver)
    N = length(names)
    for (n, name) in enumerate(names)
        field = getproperty(solver, name)
        if field !== nothing
            print(io, "\n", pad, "   ‣ " * string(name) * " ≔ ")
            show(IOContext(io, :pad => "   "), field)
        end
    end
end

"""
    BackwardEuler(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver
    ImplicitEuler(args...; kwargs...) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 1st-order backward Euler method.
"""
function BackwardEuler(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        1 1;
        1 1;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end
@doc (@doc BackwardEuler) ImplicitEuler(args...; kwargs...) = BackwardEuler(args...; kwargs...)

"""
    ImplicitMidpoint(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 2nd-order implicit midpoint method.
"""
function ImplicitMidpoint(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        1/2 1/2;
         2   1 ;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end

"""
    CrankNicolson(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 2nd-order Crank-Nicolson method.
"""
function CrankNicolson(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        0   0   0 ;
        1  1/2 1/2;
        2  1/2 1/2;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end

"""
    SDIRK3(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 3rd-order SDIRK method.
"""
function SDIRK3(; h = 0.0, ϵ = 1e-3, K = 10)
    γ = 1/2 + √3/6
    tableau = ButcherTableau(typeof(h)[
         γ    γ   0 ;
        1-γ 1-2γ  γ ;
         3   1/2 1/2;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end

"""
    GaussLegendre4(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver
    GL4(args...; kwargs...) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 4th-order Gauss-Legendre method.
"""
function GaussLegendre4(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        1/2-√3/6   1/4    1/4-√3/6;
        1/2+√3/6 1/4+√3/6   1/4   ;
           4       1/2      1/2   ;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end
@doc (@doc GaussLegendre4) GL4(args...; kwargs...) = GaussLegendre4(args...; kwargs...)

"""
    GaussLegendre6(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver
    GL6(args...; kwargs...) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 6th-order Gauss–Legendre method.
"""
function GaussLegendre6(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        1/2-√15/10     5/36    2/9-√15/15 5/36-√15/30;
           1/2     5/36+√15/24    2/9     5/36-√15/24;
        1/2+√15/10 5/36+√15/30 2/9+√15/15     5/36   ;
            6          5/18       4/9         5/18   ;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end
@doc (@doc GaussLegendre6) GL6(args...; kwargs...) = GaussLegendre6(args...; kwargs...)

"""
    LobattoIIIA4(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 4th-order Lobatto IIIA method.
"""
function LobattoIIIA4(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
         0   0    0    0  ;
        1/2 5/24 1/3 -1/24;
         1  1/6  2/3  1/6 ;
         4  1/6  2/3  1/6 ;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end

"""
    LobattoIIIB2(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 2nd-order Lobatto IIIB method.
"""
function LobattoIIIB2(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        1/2 1/2  0 ;
        1/2 1/2  0 ;
         2  1/2 1/2;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end

"""
    LobattoIIIB4(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 4th-order Lobatto IIIB method.
"""
function LobattoIIIB4(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
         0  1/6 -1/6  0 ;
        1/2 1/6  1/3  0 ;
         1  1/6  5/6  0 ;
         4  1/6  2/3 1/6;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end

"""
    LobattoIIIC2(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 2nd-order Lobatto IIIC method.
"""
function LobattoIIIC2(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        0  1/2 -1/2;
        1  1/2  1/2;
        2  1/2  1/2;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end

"""
    LobattoIIIC4(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 4th-order Lobatto IIIC method.
"""
function LobattoIIIC4(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
         0   1/6 -1/3   1/6 ;
        1/2  1/6  5/12 -1/12;
         1   1/6  2/3   1/6 ;
         4   1/6  2/3   1/6 ;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end

"""
    RadauIA3(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 3rd-order Radau IA method.
"""
function RadauIA3(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
         0  1/4 -1/4 ;
        2/3 1/4  5/12;
         3  1/4  3/4 ;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end

"""
    RadauIA5(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 5th-order Radau IA method.
"""
function RadauIA5(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
            0     1/9   -1/18-√6/18     -1/18+√6/18  ;
        3/5-√6/10 1/9  11/45+7*√6/360 11/45-43*√6/360;
        3/5+√6/10 1/9 11/45+43*√6/360  11/45-7*√6/360;
            5     1/9    4/9+√6/36       4/9-√6/36   ;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end

"""
    RadauIIA3(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 3rd-order Radau IIA method.
"""
function RadauIIA3(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        1/3  5/12 -1/12;
         1   3/4   1/4 ;
         3   3/4   1/4 ;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end

"""
    RadauIIA5(; h = 0.0, ϵ = 1e-3, K = 10) :: ImplicitRungeKuttaSolver

returns an `ImplicitRungeKuttaSolver` for the 5th-order Radau IIA method.
"""
function RadauIIA5(; h = 0.0, ϵ = 1e-3, K = 10)
    tableau = ButcherTableau(typeof(h)[
        2/5-√6/10   11/45-7*√6/360   37/225-169*√6/1800 -2/225+√6/75;
        2/5+√6/10 37/225+169*√6/1800   11/45+7*√6/360   -2/225-√6/75;
           1           4/9-√6/36          4/9+√6/36          1/9    ;
           5           4/9-√6/36          4/9+√6/36          1/9    ;
    ])
    newton = NewtonParameters(ϵ=ϵ, K=K)
    return IRK(tableau, h, newton)
end

"""
    ImplicitRungeKuttaCache{n_T, m_T, v_T, k_T, Δk_T, J_T} <: RungeKuttaCache

returns a constructor containing the temp objects needed of an `ImplicitRungeKuttaSolver`.

---

    ImplicitRungeKuttaCache(n, m, v, k, Δk, J)

returns an `ImplicitRungeKuttaCache` with:
- `n  :: Integer`                                               : step counter.
- `m  :: Integer`                                               : adaptive correction counter.
- `v  :: AbstractVector{Union{Number, AbstractVector{Number}}}` : temp for `solution.u[n]`.
- `k  :: AbstractVector{Union{Number, AbstractVector{Number}}}` : stages.
- `Δk :: AbstractVector{Union{Number, AbstractVector{Number}}}` : stages' correction.
- `J  :: AbstractMatrix{Union{Number, AbstractMatrix{Number}}}` : Jacobian of RHS function.

---

    ImplicitRungeKuttaCache(problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver)

returns an `ImplicitRungeKuttaCache` for an `ImplicitRungeKuttaSolver` given an `InitialValueProblem`.
"""
mutable struct ImplicitRungeKuttaCache{n_T, m_T, v_T, k_T, Δk_T, J_T} <: RungeKuttaCache
    n::n_T
    m::m_T
    v::v_T
    k::k_T
    Δk::Δk_T
    J::J_T
end

function ImplicitRungeKuttaCache(problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver)
    @↓ u0 = problem
    @↓ s = solver.tableau
    n = m = 1
    v = similar(u0)
    u0_T = eltype(u0)
    L = length(u0)
    k = Vector{u0_T}(undef, s, L)
    Δk = Vector{u0_T}(undef, s, L)
    J = Matrix{u0_T}(undef, L, L)
    return ImplicitRungeKuttaCache(n, m, v, k, Δk, J)
end

"""
    step!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver, cache::ImplicitRungeKuttaCache)

computes a step of the `RungeKuttaSolution` of an `InitialValueProblem` using an `ImplicitRungeKuttaSolver`.
"""
function step!(solution::RungeKuttaSolution, problem::InitialValueProblem, solver::ImplicitRungeKuttaSolver, cache::ImplicitRungeKuttaCache)
    @↓ n, v, Δk, J = cache
    @↓ u, t = solution
    k = solution.k isa Nothing ? cache.k : solution.k[n]
    @↓ f!, Df! = problem.rhs
    @↓ s, A, c, b = solver.tableau
    @↓ ϵ, K = solver.newton
    @↓ h = solver.stepsize
    v = u[n+1] # avoid allocs
    # @← J = Df(v, u[n], t[n])
    Df!(J, v, u[n], t[n])
    Z = factorize(I - h * kron(A, J))
    # Compute stages
    zero!(k)
    for l = 1:K
        for i = 1:s
            zero!(v)
            for j = 1:s
                @. v += A[i,j] * k[j]
            end
            @. v = u[n] + h * v
            # @← Δk[i] = f(v, t[n] + h * c[i])
            f!(Δk[i], v, t[n] + h * c[i])
            @. Δk[i] -= k[i]
        end
        # temporary workaround for ldiv!(Z, Δk)
        Δk_ = vcat(Δk...)
        ldiv!(Z, Δk_)
        if norm(Δk_) < ϵ * norm(k)
            break
        end
        # temporary workaround for k .+= Δk
        for i in eachindex(k)
            L = length(k[i])
            k[i] .+= Δk_[(i-1)*L+1:i*L]
        end
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
