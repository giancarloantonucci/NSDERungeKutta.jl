# NSDERungeKutta/src/phi.jl

@doc raw"""
    phifunctions(z, K::Integer; d::Integer=13) :: Vector

evaluates the exponential-integrator functions
```math
\varphi_k(z) = \frac{1}{(k-1)!} \int_0^1 e^{z(1-\theta)} \theta^{k-1} \, \mathrm{d}\theta,
\qquad \varphi_k(0) = \frac{1}{k!},
```
for ``k = 1, \dots, K``, returning them as a `Vector` in ascending order. The
algorithm is the diagonal ``(d,d)``-Padé approximation with scaling and
squaring from the EXPINT package (Berland, Skaflestad & Wright, ACM TOMS 33(1),
2007), using W. Wright's renormalised coefficient recurrences and squaring
formulae. It is a faithful port of EXPINT's `phipade.m` with these deliberate
deviations:
- No persistent global cache: in EXPINT, `phipade` memoises across calls via
  `persistent` state. Here the SOLVER CACHE owns the computed operators (they
  are evaluated once per solve at ``z = h L``), which is deterministic and
  thread-safe — global memoisation would be a data race under threaded
  Parareal.
- Polynomials are evaluated by an even/odd-split Horner rule rather than the
  Golub–Van Loan partitioned scheme of `mat_pol`; at the degrees used here
  (``d \le 13``) the two coincide in operation count and differ only by
  roundoff-neutral reassociation.
- ``e^z`` is recovered from the identity ``z\varphi_1(z) + I`` instead of a
  separate matrix exponential, consistent with the squaring stage.
- The default Padé degree is ``d = 13`` rather than EXPINT's 7. EXPINT's
  default trades accuracy near the scaling threshold (relative errors up to
  ``\sim 10^{-10}`` at ``\lVert z \rVert_\infty \approx 4``) for per-call
  speed, and its own highest-order scheme (`hochost4`) overrides it to 13.
  Here the coefficients are evaluated ONCE per solve, so the speed argument
  vanishes and ``d = 13`` holds the boundary error at roundoff
  (``\lesssim 3 \times 10^{-14}`` in testing).

# Arguments
- `z` : evaluation point, one of `Number`, `Diagonal` (evaluated elementwise —
  the fast path for spectral discretisations), dense `AbstractMatrix`, or
  sparse. Sparse DIAGONAL matrices are converted to `Diagonal`; general sparse
  matrices are densified, since ``\varphi_k`` of a sparse matrix is dense
  anyway. Krylov ``\varphi``-actions for genuinely large sparse operators are
  out of scope of this direct method.
- `K` : highest ``\varphi`` index required.
- `d` : degree of the diagonal Padé approximant (default 13; see above).
"""
function phifunctions end

phifunctions(z, K::Integer; d::Integer=13) = expphifunctions(z, K; d=d)[2]

@doc raw"""
    expphifunctions(z, K::Integer; d::Integer=13) :: Tuple

evaluates ``e^z`` TOGETHER with ``\varphi_1(z), \dots, \varphi_K(z)``,
returning the tuple `(ez, φ)` where `φ` is the `Vector` that
[`phifunctions`](@ref) returns. The exponential is the one carried through
the scaling-and-squaring chain, NOT a reconstruction from the identity
``e^z = z\varphi_1(z) + I``. The distinction matters under strong damping:
the reconstruction is only ABSOLUTELY accurate — once
``\lVert e^z \rVert`` falls below the unit roundoff of
``\lVert z\varphi_1(z) \rVert`` (real spectrum below about ``-37``) it
returns 0 or roundoff junk with UNBOUNDED relative error — whereas the
squaring chain starts from ``z_s\varphi_1(z_s) + I`` at the scaled argument,
where ``\lVert e^{z_s} \rVert \ge e^{-4}`` is far above roundoff, and
repeated squaring preserves relative accuracy down to the underflow
threshold (``e^{-700}`` to ``\sim 10^{-13}`` relative). Solver coefficients
built from this function damp stiff modes to their true ``e^{h\lambda}``
instead of a ``10^{-16}`` roundoff floor.

Dispatch mirrors [`phifunctions`](@ref): `Number`, `Diagonal` (elementwise),
dense `AbstractMatrix`, and sparse (diagonal sparse → `Diagonal`; general
sparse densified).
"""
function expphifunctions end

expphifunctions(z::Number, K::Integer; d::Integer=13) = _expphicore(z, K, d)
expphifunctions(z::AbstractMatrix{<:Number}, K::Integer; d::Integer=13) = _expphicore(z, K, d)

function expphifunctions(z::Diagonal, K::Integer; d::Integer=13)
    per = [_expphicore(λ, K, d) for λ in z.diag]
    ez = Diagonal([per[i][1] for i in eachindex(per)])
    φ = [Diagonal([per[i][2][k] for i in eachindex(per)]) for k = 1:K]
    return ez, φ
end

function expphifunctions(z::SparseArrays.AbstractSparseMatrixCSC, K::Integer; d::Integer=13)
    if isdiag(z)
        return expphifunctions(Diagonal(Vector(diag(z))), K; d=d)
    end
    return _expphicore(Matrix(z), K, d) # φₖ of a sparse matrix is dense; see docstring
end

# ------------------------------------------------------------------- core ---

_infnorm(z::Number) = abs(z)
_infnorm(z::AbstractMatrix) = opnorm(z, Inf)

function _expphicore(z, K::Integer, d::Integer)
    K ≥ 1 || throw(ArgumentError("`phifunctions` needs `K ≥ 1`."))
    Id = one(z)

    # Scaling: EXPINT's scaled_arg, s = max(0, nextpow2(‖z‖∞ / 4)).
    nrm = _infnorm(z)
    s = (nrm ≤ 4 || !isfinite(nrm)) ? 0 : ceil(Int, log2(nrm / 4))
    zs = z / 2.0^s

    # (d,d)-Padé at the scaled argument.
    N, D = padecoefficients(d, K)
    Z² = zs * zs
    φ = Vector{typeof(Z²)}(undef, K)
    for k = 1:K
        φ[k] = _evalpoly_evenodd(zs, Z², Id, D[k]) \ _evalpoly_evenodd(zs, Z², Id, N[k])
    end

    # Undo the scaling: Wright's φ-squaring recurrences (EXPINT's square_pade).
    Ez = zs * φ[1] + Id
    for _ = 1:s
        φ² = similar(φ)
        φ²[1] = (Ez + Id) * φ[1] / 2
        for k = 2:K
            v = φ[fld(k, 2)] * φ[cld(k, 2)]
            c = 2.0
            for j = k:-1:(fld(k, 2) + 1 + (isodd(k) ? 1 : 0))
                v += c * φ[j]
                c /= (k + 1 - j)
            end
            if isodd(k)
                v += φ[fld(k, 2) + 1] / factorial(fld(k, 2))
            end
            φ²[k] = v / 2.0^k
        end
        φ = φ²
        Ez = Ez * Ez
    end
    return Ez, φ
end

# Σⱼ coeff[j+1] zʲ via the even/odd split p(z) = pₑ(z²) + z⋅pₒ(z²), each part
# by Horner in z². Generic over Number and AbstractMatrix.
function _evalpoly_evenodd(z, Z², Id, coeff::AbstractVector{<:Real})
    E = _horner(Z², Id, coeff[1:2:end])
    O = _horner(Z², Id, coeff[2:2:end])
    return E + z * O
end

function _horner(Z², Id, c::AbstractVector{<:Real})
    p = c[end] * Id
    for j = length(c)-1:-1:1
        p = p * Z² + c[j] * Id
    end
    return p
end

@doc raw"""
    padecoefficients(d::Integer, K::Integer)

returns the numerator and denominator coefficient vectors (each of length
``d + 1``, constant term first) of the renormalised diagonal ``(d,d)``-Padé
approximants to ``\varphi_1, \dots, \varphi_K``, via W. Wright's recurrences
as implemented in EXPINT's `pade_cof`. The common normalisation cancels in the
division; the value at 0 is ``\varphi_k(0) = 1/k!`` by construction. All
arithmetic is floating-point from the outset: the leading constant
``(2d+1)!/d!`` overflows `Int64` already at ``d = 13``.
"""
function padecoefficients(d::Integer, K::Integer)
    n1 = prod(Float64(d + 1):Float64(2d + 1)) # (2d + 1)! / d!
    d1 = n1
    N = Vector{Vector{Float64}}(undef, K)
    D = Vector{Vector{Float64}}(undef, K)
    for ℓ = 1:K
        # Numerator: nᵢ = Σⱼ₌₀..ᵢ aᵢⱼ with a_{i0} = n1 ⋅ ∏ₘ₌₁..ᵢ 1/(ℓ+m) and
        # ratio recurrence aᵢⱼ = a_{i,j-1} ⋅ (-(d+1-j)(ℓ+1+i-j)) / ((2d+ℓ+1-j)j).
        n = Vector{Float64}(undef, d + 1)
        for i = 0:d
            a = n1
            for m = 1:i
                a /= (ℓ + m)
            end
            acc = a
            for j = 1:i
                a *= -(d + 1 - j) * (ℓ + 1 + i - j) / ((2d + ℓ + 1 - j) * j)
                acc += a
            end
            n[i+1] = acc
        end
        # Denominator: dᵢ = d1 ⋅ ∏ₘ₌₁..ᵢ (-(d+1-m)) / (m(2d+ℓ+1-m)).
        dd = Vector{Float64}(undef, d + 1)
        a = d1
        dd[1] = a
        for i = 1:d
            a *= -(d + 1 - i) / (i * (2d + ℓ + 1 - i))
            dd[i+1] = a
        end
        N[ℓ] = n
        D[ℓ] = dd
        n1 *= (2d + ℓ + 1) / (ℓ + 1)
        d1 *= (2d + ℓ + 1)
    end
    return N, D
end
