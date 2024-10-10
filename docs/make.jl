using Migmatites
using Documenter

DocMeta.setdocmeta!(Migmatites, :DocTestSetup, :(using Migmatites); recursive=true)

makedocs(;
    modules=[Migmatites],
    authors="scdyer <scdyer@uwaterloo.ca> and contributors",
    sitename="Migmatites.jl",
    format=Documenter.HTML(;
        canonical="https://sc-dyer.github.io/Migmatites.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/sc-dyer/Migmatites.jl",
    devbranch="main",
)
