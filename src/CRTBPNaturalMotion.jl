module CRTBPNaturalMotion

using StaticArrays
using OrdinaryDiffEq: ContinuousCallback, Vern9, solve

include("type_flags.jl")
include("utils.jl")

include("equations_of_motion.jl")
include("integration.jl")

export Time, ArcLength
export propagate_return_all_states, propagate_return_final_state, propagate_return_final_stm
export propagate_return_arc_length, propagate_return_time

end
