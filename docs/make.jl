using Documenter
using NSDERungeKutta

PAGES = [
    "Home" => "index.md",
    "Solvers" => "solvers.md",
    "Examples" => "examples.md",
    "API" => "api.md"
]

makedocs(;
    sitename = "NSDERungeKutta.jl",
    format = Documenter.HTML(),
    modules = [NSDERungeKutta],
    pages = PAGES,
    checkdocs = :exports, # every export must carry a docstring, or the build fails
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>"
)

deploydocs(;
    repo = "github.com/giancarloantonucci/NSDERungeKutta.jl.git"
)
