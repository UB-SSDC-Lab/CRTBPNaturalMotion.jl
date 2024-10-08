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
        "Examples" => [
            "Computing Halo Orbits" => "periodic_orbit_example.md",
            "Querying the JPL API" => "jpl_api_example.md",
        ],
        "Reference" => [
            "Public API" => "lib/public.md",
            "Internals" => "lib/internal/temp.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/UB-SSDC-Lab/CRTBPNaturalMotion.jl",
    devbranch="main",
)
