using GravityModels
using Documenter

DocMeta.setdocmeta!(GravityModels, :DocTestSetup, :(using GravityModels); recursive=true)

makedocs(;
    modules=[GravityModels],
    authors="F. Capolupo",
    sitename="GravityModels.jl",
    format=Documenter.HTML(;
        canonical="https://FraCpl.github.io/GravityModels.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/FraCpl/GravityModels.jl",
    devbranch="master",
)
