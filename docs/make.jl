using Documenter
using RungeKutta

PAGES = [
    "Home" => "index.md",
    "Examples" => "examples.md",
]

makedocs(;
    sitename = "RungeKutta.jl",
    format = Documenter.HTML(),
    modules = [RungeKutta],
    pages = PAGES,
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>"
)

deploydocs(;
    repo = "https://github.com/giancarloantonucci/RungeKutta.jl"
)
