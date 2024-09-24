using CRTBPNaturalMotion
using Documenter

DocMeta.setdocmeta!(CRTBPNaturalMotion, :DocTestSetup, :(using CRTBPNaturalMotion); recursive=true)

makedocs(;
    modules=[CRTBPNaturalMotion],
    authors="Grant Hecht",
    sitename="CRTBPNaturalMotion.jl",
    format=Documenter.HTML(;
        canonical="https://UB-SSDC-Lab.github.io/CRTBPNaturalMotion.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/UB-SSDC-Lab/CRTBPNaturalMotion.jl",
    devbranch="main",
)
