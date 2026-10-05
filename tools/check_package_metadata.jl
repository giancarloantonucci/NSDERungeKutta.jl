# Run immediately before Pkg.develop in CI. This checks the actual checkout,
# not the registered release or the active docs environment.
using TOML

root = normpath(joinpath(@__DIR__, ".."))
project = joinpath(root, "Project.toml")
println("Package root: ", root)
println("Working directory: ", pwd())
println("Active environment: ", Base.active_project())
println("CI checkout SHA: ", get(ENV, "GITHUB_SHA", "not supplied"))
isfile(project) || error("Package Project.toml is missing: $project")
println("Project.toml contents:\n", read(project, String))
metadata = TOML.parsefile(project)
get(metadata, "name", nothing) == "NSDERungeKutta" ||
    error("Expected top-level name = \"NSDERungeKutta\" before any TOML table.")
get(metadata, "uuid", nothing) == "3978d399-9848-4322-9a48-ca26b732c7dc" ||
    error("Project.toml does not have NSDERungeKutta's registered UUID.")
haskey(metadata, "version") || error("Package Project.toml has no version.")
VersionNumber(metadata["version"])
println("Package metadata checks passed.")
