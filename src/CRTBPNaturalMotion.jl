module CRTBPNaturalMotion

using Reexport
@reexport using StaticArrays
@reexport using OrdinaryDiffEq

using FastChebInterp
using Format
using LinearAlgebra
using LinearSolve
using NonlinearSolve

# For JPL API
import HTTP
import JSON3
import StructTypes

include("type_flags.jl")
include("utils.jl")
include("interpolation.jl")

include("equations_of_motion.jl")
include("integration.jl")
include("periodic_orbit_shooting_functions.jl")
include("periodic_orbit_correction.jl")
include("periodic_orbits.jl")
include("manifolds.jl")

# JPL API
include(joinpath("jpl_catalog", "query_construction.jl"))
include(joinpath("jpl_catalog", "query_structs.jl"))
include(joinpath("jpl_catalog", "query_request.jl"))
include(joinpath("jpl_catalog", "get_jpl_orbits.jl"))

# Provides easy definition of sample orbits (not exported)
include("scenario_orbits.jl")

# Type flags and utils
export Time, ArcLength
export jacobi_integral

# Interpolation functions
export value, jacobian

# Propagation functions
export propagate_return_all_states, propagate_return_final_state, propagate_return_final_stm
export propagate_return_arc_length, propagate_return_time

# Periodic orbits
export correct_typeA_initial_conditions
export TypeAPeriodicOrbit, GeneralPeriodicOrbit, get_full_orbit
export generate_periodic_orbit_cheb_interpolation
export generate_periodic_orbit_cheb_approximation

# Invariant manifolds
export InvariantManifold
export get_stable_manifold_trajectory
export get_unstable_manifold_trajectory

# JPL API
export get_jpl_orbits
export minimum_period_orbit, maximum_period_orbit

end
