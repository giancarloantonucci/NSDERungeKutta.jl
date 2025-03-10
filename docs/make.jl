using Documenter
using DocumenterInterLinks
using NSDEBase, NSDERungeKutta

PAGES = [
    "Home" => "index.md",
    "Examples" => "examples.md",
]

links = InterLinks(
    "NSDEBase" => (
        "https://giancarloantonucci.github.io/NSDEBase.jl/dev/",
        "https://giancarloantonucci.github.io/NSDEBase.jl/dev/objects.inv"
    )
)

makedocs(;
    sitename = "NSDERungeKutta.jl",
    format = Documenter.HTML(),
    modules = [NSDERungeKutta],
    pages = PAGES,
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>",
    plugins = [links],
)

deploydocs(;
    repo = "https://github.com/giancarloantonucci/NSDERungeKutta.jl"
)
