module CRTBPNaturalMotion

using Reexport
@reexport using StaticArrays
@reexport using OrdinaryDiffEq

using FastChebInterp
using LinearAlgebra
using LinearSolve
using NonlinearSolve

include("type_flags.jl")
include("utils.jl")
include("interpolation.jl")

include("equations_of_motion.jl")
include("integration.jl")
include("periodic_orbit_shooting_functions.jl")
include("periodic_orbit_correction.jl")
include("periodic_orbits.jl")
include("manifolds.jl")

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
export TypeAPeriodicOrbit, get_full_orbit
export generate_periodic_orbit_cheb_interpolation
export generate_periodic_orbit_cheb_approximation

# Invariant manifolds
export InvariantManifold
export get_stable_manifold_trajectory
export get_unstable_manifold_trajectory

end
