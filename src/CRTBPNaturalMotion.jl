module CRTBPNaturalMotion

using Reexport
@reexport using StaticArrays
@reexport using OrdinaryDiffEq
using NonlinearSolve
using LinearSolve

include("type_flags.jl")
include("utils.jl")

include("equations_of_motion.jl")
include("integration.jl")
include("periodic_orbit_shooting_functions.jl")
include("periodic_orbit_correction.jl")
include("PeriodicOrbitStructs.jl")

# Provides easy definition of sample orbits (not exported)
include("scenario_orbits.jl")

# Type flags and utils
export Time, ArcLength
export jacobi_integral

# Propagation functions
export propagate_return_all_states, propagate_return_final_state, propagate_return_final_stm
export propagate_return_arc_length, propagate_return_time

# Periodic orbits
export correct_typeA_initial_conditions
export TypeAPeriodicOrbit, get_full_orbit

end
