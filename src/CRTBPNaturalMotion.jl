module CRTBPNaturalMotion

using Reexport
@reexport using StaticArrays
@reexport using OrdinaryDiffEq

using JLD2: load, save
using FastChebInterp: chebregression, chebinterp, chebpoints, chebjacobian, FastChebInterp
using Format: printfmt
using LinearAlgebra: eigen, norm
using NonlinearSolve: NonlinearFunction, NonlinearProblem, SimpleTrustRegion

# For AD
import DifferentiationInterface as AD
using Enzyme: Enzyme

# For JPL API
using HTTP: HTTP
using JSON3: JSON3
using StructTypes: StructTypes

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

import PrecompileTools
PrecompileTools.@compile_workload begin
    # Define some constants to use in precompilation
    G = 6.673e-20
    m_1 = 5.9742e24         # Mass of the Earth in kg
    m_2 = 7.3483e22         # Mass of the Moon in kg
    r_12 = 384400.0         # Mean distance between the Earth and Moon
    mu = m_2 / (m_1 + m_2)   # The mass ratio for the system
	DU = r_12
	TU = 1.0 / sqrt((G * (m_1 + m_2)) / DU^3)

    # ==== Common propagation functions
    x0 = SA[0.823, 0.0, -0.022, 0.0, 0.134, 0.0]
    ts = (0.0, 1.0/TU)
    ls = (0.0, 1.0/DU)
    propagate_return_all_states(x0, ts, mu)
    propagate_return_all_states(x0, ts, mu, ArcLength)
    propagate_return_final_state(x0, ts, mu)
    propagate_return_final_state(x0, ts, mu, ArcLength)

    # ==== Type A periodic orbits
    rx_guess = 0.8233851820
    rz_guess = -0.0222775563
    vy_guess = 0.1341841703
    C_guess = jacobi_integral(SA[rx_guess, 0.0, rz_guess, 0.0, vy_guess, 0.0], mu)
    typeA = TypeAPeriodicOrbit(
        rx_guess, rz_guess, vy_guess, mu, TU, DU;
        N_cross = 1,
        constraint = (:jacobi_integral, C_guess - 0.001),
    )

    # ==== JPL Query/GeneralPeriodicOrbit
    gen = GeneralPeriodicOrbit(
        SA[
            0.8278817186555927, 1.0963203778218452e-27, 0.09920318722755168,
            1.8521961003574663e-15, 0.2146176871613027, 6.056332850673628e-15,
        ],
        3.10413141229191, 300.038529354993, 2.785539896020308, 0.6175420021307126,
        0.01215058560962404, 382981.289129055, 389703.264829278,
    )

    # ==== Orbit Trajectories and Approximations
    for orbit in (typeA, gen)
        get_full_orbit(orbit)

        orbit_interp = generate_periodic_orbit_cheb_interpolation(
            orbit, 10, Time,
        )
        orbit_interp = generate_periodic_orbit_cheb_interpolation(
            orbit, 10, ArcLength,
        )
        value(orbit_interp, 0.5)
        derivative(orbit_interp, 0.5)
        second_derivative(orbit_interp, 0.5)

        orbit_lsqf = generate_periodic_orbit_cheb_approximation(
            orbit, 50, 10, Time,
        )
        orbit_lsqf = generate_periodic_orbit_cheb_approximation(
            orbit, 50, 10, ArcLength,
        )
        value(orbit_lsqf, 0.5)
        derivative(orbit_lsqf, 0.5)
        second_derivative(orbit_lsqf, 0.5)

        man = InvariantManifold(orbit, 50.0 / DU)
        man_start = 0.0
        man_length = 10.0/DU

        man_cs_interp = generate_stable_manifold_cross_section_cheb_interpolant(
            man, true, 10, ArcLength, 0.5, ArcLength;
            start_cond = man_start,
            manifold_length = man_length,
        )
        value(man_cs_interp, 0.5)
        derivative(man_cs_interp, 0.5)
        second_derivative(man_cs_interp, 0.5)

        man_cs_lsqf = generate_stable_manifold_cross_section_cheb_approximation(
            man, true, 50, 10, ArcLength, 0.5, ArcLength;
            start_cond = man_start,
            manifold_length = man_length,
        )
        value(man_cs_lsqf, 0.5)
        derivative(man_cs_lsqf, 0.5)
        second_derivative(man_cs_lsqf, 0.5)

        man_interp = generate_stable_manifold_cheb_interpolant(
            man, true, 10, ArcLength, 10, ArcLength;
            start_cond = man_start,
            manifold_length = man_length,
        )
        value(man_interp, 0.5, 0.5)
        jacobian(man_interp, 0.5, 0.5)
        hessian(man_interp, 0.5, 0.5)

        man_lsqf = generate_stable_manifold_cheb_approximation(
            man, true, 10, 5, ArcLength, 10, 5, ArcLength;
            start_cond = man_start,
            manifold_length = man_length,
        )
        value(man_lsqf, 0.5, 0.5)
        jacobian(man_lsqf, 0.5, 0.5)
        hessian(man_lsqf, 0.5, 0.5)
    end
end

# Type flags and utils
export Time, ArcLength
export jacobi_integral

# Interpolation functions
export value, derivative, second_derivative, jacobian, hessian

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
export generate_stable_manifold_cross_section_cheb_interpolant
export generate_stable_manifold_cross_section_cheb_approximation
export generate_stable_manifold_cheb_interpolant
export generate_stable_manifold_cheb_approximation
export save_interp, load_interp

# JPL API
export get_jpl_orbits
export minimum_period_orbit, maximum_period_orbit

end
