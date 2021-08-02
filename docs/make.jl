using Documenter
using RungeKutta

PAGES = ["Home" => "index.md"]

makedocs(;
    sitename = "RungeKutta",
    format = Documenter.HTML(),
    modules = [RungeKutta],
    pages = PAGES,
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>"
)

deploydocs(;
    repo = "https://github.com/antonuccig/RungeKutta.jl"
)
