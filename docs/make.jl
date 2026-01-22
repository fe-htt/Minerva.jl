using Minerva
using Documenter

DocMeta.setdocmeta!(Minerva, :DocTestSetup, :(using Minerva); recursive=true)

makedocs(;
    modules=[Minerva],
    authors="Federico Hattab",
    sitename="Minerva.jl",
    format=Documenter.HTML(;
        canonical="https://fe-htt.github.io/Minerva.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fe-htt/Minerva.jl",
    devbranch="main",
)
