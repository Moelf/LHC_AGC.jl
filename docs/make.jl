using Documenter
using Pluto: Configuration.CompilerOptions
using PlutoStaticHTML

notebooks = [
    "Notebook",
]

include("build.jl")

build()
md_files = markdown_files()
T = [t => f for (t, f) in zip(notebooks, md_files)]

makedocs(;
    modules=Module[],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    pages=[
        "Notebooks" => T,
    ],
    repo="https://github.com/Moelf/LHC_AGC.jl/blob/{commit}{path}#L{line}",
    sitename="LHC Analysis Grand Challengt",
    authors="Jerry Ling",
    assets=String[],
)

deploydocs(;
    repo="github.com/Moelf/LHC_AGC.jl",
)
