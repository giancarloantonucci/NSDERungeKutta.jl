# NSDERungeKutta/test/oracle.jl
#
# Mathematical verification oracle for the exprk solvers — include from
# runtests.jl AFTER exprk.jl. Adapted from the standalone ExpRKOracle.jl:
# @test instead of the script's gate() (so failures FAIL Pkg.test), testset
# scoping instead of globals, helper names deconflicted. Costs ~6-8 s.
#
# References: [KT05] Kassam & Trefethen SISC 26(4) 2005 (the naive-formula
# instability near zero eigenvalues, contour-integral fix); [TW14] Trefethen
# & Weideman SIREV 56(3) 2014 (exponentially convergent trapezoidal rule for
# Dunford-Taylor integrals); [BSW07] Berland-Skaflestad-Wright ACM TOMS 33(1)
# 2007 (phipade, the stable Padé route this package ships); [HO10] Hochbruck
# & Ostermann Acta Numerica 19 2010 (the ψ order conditions).
#
# Two evaluation routes INDEPENDENT of the shipped Padé engine — 256-bit
# direct formulas for scalars, [TW14] contour integrals for matrices — plus
# the exact [HO10] operator identities of every tableau, plus the [KT05]
# Kuramoto-Sivashinsky benchmark whose symbol carries an exact zero
# eigenvalue. The identity gate has a two-for-two record as a bug detector:
# ψ₁ convicted EXPINT's shipped ETD2RK weights, and the relative U/V
# comparison caught the zφ₁ + I exponential reconstruction (absolutely
# accurate, relatively meaningless under strong damping) that
# expphifunctions now replaces.

using Test, LinearAlgebra, Printf

# ---- oracle 1: 256-bit direct formulas (unusable in doubles — that is the
# point: at 256 bits the cancellation costs ~77 digits of headroom and
# leaves ~60 correct, so it is ground truth, not an algorithm) -------------
function oracle_φbig(z::Number, K::Integer)
    setprecision(BigFloat, 256) do
        zb = Complex{BigFloat}(z)
        if abs(zb) < big"1e-40"
            return [ComplexF64(1 / factorial(big(k))) for k = 1:K]
        end
        et = exp(zb)
        part = one(zb)
        tk = zb
        out = Vector{ComplexF64}(undef, K)
        for k = 1:K
            out[k] = ComplexF64((et - part) / tk)
            part += tk / factorial(big(k))
            tk *= zb
        end
        return out
    end
end

# ---- oracle 2: [TW14] trapezoidal contour integral; pointwise φ on the
# contour by a stable hybrid (series for |t| ≤ 1, direct beyond) -----------
function oracle_φhybrid(t::Complex, K::Integer)
    out = Vector{ComplexF64}(undef, K)
    if abs(t) ≤ 1
        for k = 1:K
            s = zero(ComplexF64)
            term = ComplexF64(1 / factorial(k))
            for j = 0:29
                s += term
                term *= t / (j + k + 1)
            end
            out[k] = s
        end
    else
        et = exp(t)
        part = one(ComplexF64)
        tk = t
        for k = 1:K
            out[k] = (et - part) / tk
            part += tk / factorial(k)
            tk *= t
        end
    end
    return out
end

function oracle_φcontour(Z::AbstractMatrix, K::Integer; N::Integer=512, margin::Real=3.0)
    ev = eigvals(Matrix{ComplexF64}(Z))
    c = (minimum(real, ev) + maximum(real, ev)) / 2 + im * (minimum(imag, ev) + maximum(imag, ev)) / 2
    r = maximum(abs, ev .- c) + margin
    n = size(Z, 1)
    Id = Matrix{ComplexF64}(I, n, n)
    acc = [zeros(ComplexF64, n, n) for _ = 1:K]
    for j = 0:N-1
        w = r * cis(2π * (j + 0.5) / N)
        t = c + w
        R = (t * Id - Z) \ Id
        f = oracle_φhybrid(t, K)
        wR = w .* R
        for k = 1:K
            acc[k] .+= f[k] .* wR
        end
    end
    return [a ./ N for a in acc]
end

oracle_relerr(a, b) = norm(a - b) / max(norm(b), eps())

# Chebyshev differentiation matrix (Trefethen, Spectral Methods in MATLAB,
# p. 54) — the non-normal Allen-Cahn operator of [KT05] §5.
function oracle_cheb(N::Integer)
    x = [cos(π * i / N) for i = 0:N]
    c = [(i == 0 || i == N ? 2.0 : 1.0) * (-1.0)^i for i = 0:N]
    D = [i == j ? 0.0 : c[i+1] / c[j+1] / (x[i+1] - x[j+1]) for i = 0:N, j = 0:N]
    D -= Diagonal(vec(sum(D, dims=2)))
    return D
end

# The [KT05] KS symbol on [0, 32π] with d Fourier modes: exact zero
# eigenvalue at κ = 0, small POSITIVE ones at the lowest modes.
function oracle_ks_spectrum(d::Integer)
    m = [j ≤ d ÷ 2 ? j - 1 : j - 1 - d for j = 1:d]
    κ = m ./ 16
    return κ, κ .^ 2 .- κ .^ 4
end

oracle_solvers = (LawsonEuler, NorsettEuler, ETD2RK, ETD3RK, ETD4RK, Lawson4, Krogstad, HochbruckOstermann4)
oracle_Zcheb = 0.25 .* (0.01 .* (oracle_cheb(20) * oracle_cheb(20))[2:end-1, 2:end-1])

@testset "oracle: φ engine vs 256-bit ground truth [KT05 poison set]" begin
    _, λ128 = oracle_ks_spectrum(128)
    poison = ComplexF64[1e-2, 1e-6, 1e-8, -1e-8, 0.0, -0.5, 3.9, -40.0, 137.0, 2.5 + 1.3im]
    append!(poison, ComplexF64.(0.25 .* λ128[1:4]))
    push!(poison, ComplexF64(0.25 * minimum(λ128)))
    worst = 0.0
    for z in poison
        φp = imag(z) == 0 ? phifunctions(real(z), 3) : phifunctions(z, 3)
        φb = oracle_φbig(z, 3)
        for k = 1:3
            worst = max(worst, abs(ComplexF64(φp[k]) - φb[k]) / abs(φb[k]))
        end
    end
    @test worst < 1e-12
    # chain exponential under strong damping (RELATIVE accuracy):
    for z in (-40.0, -200.0, -700.0)
        ez, = expphifunctions(z, 1)
        @test abs(ez - exp(z)) / exp(z) < 1e-11
    end
    # [KT05] Table 2.2 re-enactment (γ = 4φ₃ - φ₂; the naive formula returns
    # -888.17 at z = 1e-6 where the truth is 1/6):
    for z in (1e-4, 1e-6)
        φ = phifunctions(z, 3)
        γtruth = real(4 * oracle_φbig(z, 3)[3] - oracle_φbig(z, 3)[2])
        @test abs((4φ[3] - φ[2]) - γtruth) / abs(γtruth) < 1e-13
    end
end

@testset "oracle: φ engine vs [TW14] contour integral on matrices" begin
    Msym = [sin(3i + j) for i = 1:6, j = 1:6] - 3.0I
    Msym = (Msym + Msym') / 2
    for Z in (Msym, oracle_Zcheb) # dense symmetric; non-normal Chebyshev z = h⋅L
        φp = phifunctions(Z, 3)
        φc = oracle_φcontour(Z, 3; N=512)
        φc2 = oracle_φcontour(Z, 3; N=1024)
        @test maximum(oracle_relerr(φp[k], φc[k]) for k = 1:3) < 1e-11
        @test maximum(oracle_relerr(φc[k], φc2[k]) for k = 1:3) < 1e-12 # oracle self-consistency
    end
end

@testset "oracle: [HO10] operator identities of the shipped tableaus" begin
    ψ_applicable = Dict(:LawsonEuler => 0, :NorsettEuler => 1, :ETD2RK => 2, :ETD3RK => 3,
                        :ETD4RK => 3, :Lawson4 => 0, :Krogstad => 3, :HochbruckOstermann4 => 3)
    rows_are_φform = Dict(:LawsonEuler => false, :NorsettEuler => true, :ETD2RK => true, :ETD3RK => true,
                          :ETD4RK => true, :Lawson4 => false, :Krogstad => true, :HochbruckOstermann4 => true)
    _, λ128 = oracle_ks_spectrum(128)
    testargs = Any[hcat(z) for z in (1e-8, 1e-2, -40.0, 137.0)]
    push!(testargs, hcat(2.5 + 1.3im))
    push!(testargs, Diagonal(0.25 .* λ128)) # exact zero + near-zero positive + stiff
    push!(testargs, oracle_Zcheb)           # non-normal dense
    worstψ = worstrow = worstexp = 0.0
    for maker in oracle_solvers
        tab = maker(h=1.0).tableau
        nψ = ψ_applicable[tab.name]
        for z in testargs
            U, V, A, B = tab.κ(z)
            φ = phifunctions(z, 3)
            scale = max(norm(φ[1]), 1.0)
            # U/V exponential identities — RELATIVE, elementwise on spectra:
            expref(w) = w isa Diagonal ? Diagonal(exp.(w.diag)) : exp(Matrix(w))
            exprel(a, ref) = ref isa Diagonal ?
                maximum(abs.((a isa Diagonal ? a.diag : diag(Matrix(a))) .- ref.diag) ./ max.(abs.(ref.diag), floatmin())) :
                oracle_relerr(Matrix(a), ref)
            ez = expref(z)
            worstexp = max(worstexp, exprel(V, ez))
            for i = 1:tab.s
                ci = tab.c[i]
                target = ci == 0 ? (z isa Diagonal ? Diagonal(ones(size(z, 1))) : Matrix(one(z))) :
                         ci == 1 ? ez : expref(ci * z)
                worstexp = max(worstexp, exprel(U[i], target))
            end
            # ψ weight identities:
            if nψ ≥ 1
                worstψ = max(worstψ, norm(sum(B[i] for i = 1:tab.s if B[i] !== nothing) - φ[1]) / scale)
            end
            if nψ ≥ 2
                worstψ = max(worstψ, norm(sum(B[i] * tab.c[i] for i = 1:tab.s if B[i] !== nothing) - φ[2]) / max(norm(φ[2]), 1.0))
            end
            if nψ ≥ 3
                worstψ = max(worstψ, norm(sum(B[i] * tab.c[i]^2 for i = 1:tab.s if B[i] !== nothing) / 2 - φ[3]) / max(norm(φ[3]), 1.0))
            end
            # internal-row consistency (φ-form schemes) or z = 0 consistency (Lawson):
            if rows_are_φform[tab.name]
                for i = 2:tab.s
                    rowsum = sum(A[i, j] for j = 1:i-1 if A[i, j] !== nothing)
                    worstrow = max(worstrow, norm(rowsum - tab.c[i] * phifunctions(tab.c[i] * z, 1)[1]) / scale)
                end
            else
                z0 = zero(Matrix(z))
                U0, V0, A0, B0 = tab.κ(z0)
                worstψ = max(worstψ, norm(sum(B0[i] for i = 1:tab.s if B0[i] !== nothing) - one(z0)))
                for i = 2:tab.s
                    rowsum = sum(A0[i, j] for j = 1:i-1 if A0[i, j] !== nothing)
                    worstrow = max(worstrow, norm(rowsum - tab.c[i] * one(z0)))
                end
            end
        end
    end
    @test worstψ < 1e-12
    @test worstrow < 1e-12
    @test worstexp < 1e-11
end

@testset "oracle: [KT05] Kuramoto-Sivashinsky benchmark" begin
    # u_t = -u u_x - u_xx - u_xxxx on [0, 32π], Fourier space, dense DFT
    # matrices (dependency-free); exact zero and unstable near-zero
    # eigenvalues in L — the regime of [KT05] Fig. 3.5. Order against a
    # cross-family RK4 reference.
    dks = 128
    κ, λ = oracle_ks_spectrum(dks)
    W = [cis(-2π * m * j / dks) for m = 0:dks-1, j = 0:dks-1]
    Winv = W' ./ dks
    xs = [32π * j / dks for j = 0:dks-1]
    û0 = W * complex.(@. cos(xs / 16) * (1 + sin(xs / 16)))
    Lks = Diagonal(complex.(λ))
    nonstiff(û, t) = (-im / 2) .* κ .* (W * ((Winv * û) .^ 2))
    ksproblem = IVP(SRHS(Lks, RHS(nonstiff; iscomplex=true)), û0, (0.0, 1.0))
    fullrhs(û, t) = Lks * û .+ nonstiff(û, t)
    uref = Winv * solve(IVP(RHS(fullrhs; iscomplex=true), û0, (0.0, 1.0)), RK4(h=1e-3)).u[end]
    for maker in (ETD4RK, HochbruckOstermann4)
        hs = [0.1, 0.05, 0.025, 0.0125]
        errs = map(hs) do h
            norm(Winv * solve(ksproblem, maker(h=h)).u[end] - uref) / norm(uref)
        end
        lx, ly = log.(hs), log.(errs)
        order = sum((lx .- sum(lx) / 4) .* (ly .- sum(ly) / 4)) / sum(abs2, lx .- sum(lx) / 4)
        @test order > 3.5
        @test issorted(errs; rev=true)
    end
end
