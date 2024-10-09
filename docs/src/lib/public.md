# Public API Documentation

Documentation for CRTBPNaturalMotion's public interface.

## Contents
```@contents
Pages = ["public.md"]
Depth = 2:2
```

## Index 
```@index
Pages = ["public.md"]
```

## Querying JPL CRTBP Periodic Orbit API
```@autodocs
Modules = [CRTBPNaturalMotion]
Pages = [joinpath("jpl_catalog","get_jpl_orbits.jl")]
Private = false
```

## Periodic Orbit Computation
```@autodocs
Modules = [CRTBPNaturalMotion]
Pages = ["periodic_orbits.jl", "periodic_orbit_correction.jl"]
Private = false
```

## Invariant Manifold Computation
```@autodocs
Modules = [CRTBPNaturalMotion]
Pages = ["manifolds.jl"]
Private = false
```

## Type Flags
```@autodocs
Modules = [CRTBPNaturalMotion]
Pages = ["type_flags.jl"]
Private = false
```

## Propagation Functions
```@autodocs
Modules = [CRTBPNaturalMotion]
Pages = ["integration.jl"]
Private = false
```

## Utility Functions
```@autodocs
Modules = [CRTBPNaturalMotion]
Pages = ["utils.jl", "interpolation.jl"]
Private = false
```

