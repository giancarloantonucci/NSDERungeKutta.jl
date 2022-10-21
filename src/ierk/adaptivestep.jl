# function adaptivestep!(cache::ImplicitExplicitRungeKuttaCache, solution::AbstractRungeKuttaSolution, solver::ImplicitExplicitRungeKuttaSolver, adaptive::AdaptiveParameters)
#     @↓ n, m, kᴵ, kᴱ = cache
#     @↓ u, savestages = solution
#     @↓ implicit, explicit, stepsize = solver
#     @↓ s, bᴵ ← b, pᴵ ← p, dᴵ ← d, qᴵ ← q = implicit.tableau
#     @↓ bᴱ ← b, pᴱ ← p, dᴱ ← d, qᴱ ← q = explicit.tableau
#     @↓ h = stepsize
#     @↓ atol, rtol, nits = adaptive
#     zero!(v)
#     for i = 1:s
#         @. v += (bᴵ[i] - dᴵ[i]) * kᴵ[i]
#         @. v += (bᴱ[i] - dᴱ[i]) * kᴱ[i]
#     end
#     tol = atol + norm(u[n]) * rtol
#     err = hairernorm(v)
#     pwr = 1 / (min(pᴱ, qᴱ) + 1)
#     r = max(0.5, min(2.0, (0.35 * tol / err) ^ pwr))
#     h *= r
#     if err < tol || m ≥ nits
#         n += 1
#     else
#         m += 1
#     end
#     @↑ cache = n, m
#     @↑ stepsize = h
#     return solution
# end
