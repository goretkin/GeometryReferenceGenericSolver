using Documenter, GeometryReferenceGenericSolver

makedocs(
    modules = [GeometryReferenceGenericSolver],
    format = Documenter.HTML(),
    checkdocs = :exports,
    sitename = "GeometryReferenceGenericSolver.jl",
    pages = Any["index.md"]
)

deploydocs(
    repo = "github.com/goretkin/GeometryReferenceGenericSolver.jl.git",
)
