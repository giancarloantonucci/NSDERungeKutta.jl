function Base.show(io::IO, tableau::ButcherTableau)
    if get(io, :compact, false)
        print(io, "ButcherTableau")
    else
        @↓ A, b, c, s, p, d, q = tableau
        print(io, "ButcherTableau:\n",
            "  A: ", A, "\n",
            "  b: ", b, "\n",
            "  c: ", c, "\n",
            "  s: ", s, "\n",
            "  p: ", p, "\n",
            "  d: ", d, "\n",
            "  q: ", q, "\n",
        )
    end
end

function Base.show(io::IO, adaptive::AdaptiveParameters)
    if get(io, :compact, false)
        print(io, "AdaptiveParameters")
    else
        @↓ δ, ϵ, K = adaptive
        print(io,
            "AdaptiveParameters:\n",
            "  δ: ", δ, "\n",
            "  ϵ: ", ϵ, "\n",
            "  K: ", K, "\n",
        )
    end
end

function Base.show(io::IO, erk::ExplicitRungeKuttaSolver)
    if get(io, :compact, false)
        print(io, "ExplicitRungeKuttaSolver")
    else
        @↓ tableau, h, adaptive = erk
        print(io,
            "ExplicitRungeKuttaSolver:\n",
            "  tableau: ",  tableau, "\n",
            "  h: ", h, "\n",
            "  adaptive: ", typeof(adaptive), "\n",
        )
    end
end

function Base.show(io::IO, irk::ImplicitRungeKuttaSolver)
    if get(io, :compact, false)
        print(io, "ImplicitRungeKuttaSolver")
    else
        @↓ tableau, h, ϵ, K, adaptive = irk
        print(io,
            "ImplicitRungeKuttaSolver:\n",
            "  tableau: ",  typeof(tableau), "\n",
            "  h: ", h, "\n",
            "  ϵ: ", ϵ, "\n",
            "  K: ", K, "\n",
            "  adaptive: ", typeof(adaptive), "\n",
        )
    end
end

function Base.show(io::IO, solution::RungeKuttaSolution)
    if get(io, :compact, false)
        print(io, "RungeKuttaSolution")
    else
        @↓ u, t, k = solution
        print(io,
            "RungeKuttaSolution:\n",
            "  u: ", u, "\n",
            "  t: ", t, "\n",
            "  k: ", k, "\n",
        )
    end
end
