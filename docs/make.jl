using CRTBPNaturalMotion
using Documenter

#DocMeta.setdocmeta!(CRTBPNaturalMotion, :DocTestSetup, :(using CRTBPNaturalMotion); recursive=true)

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
            #"Computing Halo Orbits" => "periodic_orbit_example.md",
            #"Querying the JPL API"  => "jpl_api_example.md",
            "Interpolation"         => "interpolation_example.md",
        ],
        "Reference" => [
            "Public API"                    => "lib/public.md",
            "Querying JPL Catelog"          => "lib/public/jpl_orbits.md",
            "Periodic Orbit Computation"    => "lib/public/orbit_computation.md",
            "Stable/Unstable Manifolds"     => "lib/public/manifolds.md",
            "Propagation"                   => "lib/public/propagation.md",
            "Type Flags"                    => "lib/public/type_flags.md",
            "Utilities"                     => "lib/public/utilities.md",
        ],
        "Developer" => [
            "Internals" => "lib/internal/temp.md",
        ]
    ],
    pagesonly = true,
)

deploydocs(;
    repo="github.com/UB-SSDC-Lab/CRTBPNaturalMotion.jl",
    devbranch="main",
)
