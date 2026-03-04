using Minerva
using Minerva.PostprocessingUtils
using Documenter

DocMeta.setdocmeta!(Minerva, :DocTestSetup, :(using Minerva); recursive=true)
include("pages.jl")
 
makedocs(;
    modules=[
        Minerva
        Minerva.PostprocessingUtils
    ],
    authors="Federico Hattab",
    sitename="Minerva",
    clean = true, doctest = false, #linkcheck = true,
    warnonly = [:docs_block, :missing_docs, :cross_references],
    format=Documenter.HTML(;
        canonical="https://fe-htt.github.io/Minerva.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=pages
)

deploydocs(;
    repo="github.com/fe-htt/Minerva.jl.git",
    devbranch="main",
)
