using Documenter
using DocumenterInterLinks
using NSDEBase, NSDERungeKutta

PAGES = [
    "Home" => "index.md",
    "Examples" => "examples.md",
]

# Ensure the inventories directory exists
mkpath(joinpath(@__DIR__, "inventories"))

links = InterLinks(
    "sphinx" => "https://www.sphinx-doc.org/en/master/",
    "matplotlib" => "https://matplotlib.org/3.7.3/objects.inv",
    "Julia" => (
        "https://docs.julialang.org/en/v1/",
        joinpath(@__DIR__, "inventories", "Julia.toml")
    ),
    "Documenter" => (
        "https://documenter.juliadocs.org/stable/",
        "https://documenter.juliadocs.org/stable/objects.inv",
        joinpath(@__DIR__, "inventories", "Documenter.toml")
    ),
);

makedocs(;
    sitename = "NSDERungeKutta.jl",
    format = Documenter.HTML(),
    # modules = [NSDERungeKutta],
    pages = PAGES,
    authors = "Giancarlo A. Antonucci <giancarlo.antonucci@icloud.com>",
    plugins = links,
)

deploydocs(;
    repo = "https://github.com/giancarloantonucci/NSDERungeKutta.jl"
)
