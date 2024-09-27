
abstract type AbstractPeriodicOrbit end

"""
    ParameterizationType

A Union type for the different types of parameterizations that can be used when parameterizing
a periodic orbit or its associated invairiant stable or intable manifold.
"""
ParameterizationType = Union{Type{Time}, Type{ArcLength}}

"""
    TypeAPeriodicOrbit

A representation of a CRTBP periodic orbit of Type A.

# Fields
- `u0::SVector{3, Float64}`: The initial conditions for the periodic orbit.
- `P::Float64`: The period of the periodic orbit.
- `S::Float64`: The arc length along the periodic orbit.
- `N_cross::Int`: The number of crossings through x-z plane in half period.
- `mu::Float64`: The gravitational parameter of the CRTBP.
"""
struct TypeAPeriodicOrbit <: AbstractPeriodicOrbit
    # The initial conditions for the periodic orbit
    u0::SVector{3, Float64}

    # The period of the periodic orbit
    P::Float64

    # The arc length along the periodic orbit
    S::Float64

    # The number of crossings through x-z plane in half period
    N_cross::Int

    # The gravitational parameter
    mu::Float64
end

"""
    TypeAPeriodicOrbit(
        rx, rz, vy, mu;
        N_cross = 1,
        constraint = (:x_start_coordinate, 0.8),
        ode_solver = Vern9(),
        ode_reltol = 1e-12,
        ode_abstol = 1e-12,
        nl_solver  = SimpleTrustRegion(autodiff = nothing, nlsolve_update_rule = Val(true)),
    )
"""
function TypeAPeriodicOrbit(
    rx, rz, vy, mu;
    N_cross = 1,
    constraint = (:x_start_coordinate, 0.8),
    ode_solver = Vern9(),
    ode_reltol = 1e-12,
    ode_abstol = 1e-12,
    nl_solver  = SimpleTrustRegion(autodiff = nothing, nlsolve_update_rule = Val(true)),
)
    # Correct the initial conditions
    sol, cflag =  correct_typeA_initial_conditions(
        rx, rz, vy, mu;
        N_cross = N_cross,
        constraint = constraint,
        ode_solver = ode_solver,
        ode_reltol = ode_reltol,
        ode_abstol = ode_abstol,
        nl_solver  = nl_solver,
    )
    rx, rz, vy, P = sol

    if !cflag
        @warn "Periodic orbit correction did not converge. Consider trying again with different initial guess"
    end

    # Compute orbit arc length
    S = propagate_return_arc_length(
        SA[rx, 0.0, rz, 0.0, vy, 0.0], (0.0, P), mu;
        solver = ode_solver,
        reltol = ode_reltol,
        abstol = ode_abstol,
    )

    return TypeAPeriodicOrbit(
        SA[rx,rz,vy], P, S, N_cross, mu,
    )
end
"""
    TypeAPeriodicOrbit(
        orbit::TypeAPeriodicOrbit;
        constraint = (:x_start_coordinate, 0.8),
        ode_solver = Vern9(),
        ode_reltol = 1e-12,
        ode_abstol = 1e-12,
        nl_solver  = SimpleTrustRegion(autodiff = nothing, nlsolve_update_rule = Val(true)),
    )
"""
function TypeAPeriodicOrbit(
    orbit::TypeAPeriodicOrbit;
    constraint = (:x_start_coordinate, 0.8),
    ode_solver = Vern9(),
    ode_reltol = 1e-12,
    ode_abstol = 1e-12,
    nl_solver  = SimpleTrustRegion(autodiff = nothing, nlsolve_update_rule = Val(true)),
)
    return TypeAPeriodicOrbit(
        orbit.u0[1], orbit.u0[2], orbit.u0[3], orbit.mu;
        N_cross = orbit.N_cross,
        constraint = constraint,
        ode_solver = ode_solver,
        ode_reltol = ode_reltol,
        ode_abstol = ode_abstol,
        nl_solver  = nl_solver,
    )
end

"""
    get_full_initial_state(orbit::TypeAPeriodicOrbit)

Get the full initial reference state of the periodic `orbit`.

# Arguments
- `orbit::TypeAPeriodicOrbit`: The periodic orbit.

# Returns
- `x0::SVector{6,Float64}`: The full initial state of the periodic orbit.
"""
get_full_initial_state(orbit::TypeAPeriodicOrbit) = SA[orbit.u0[1], 0.0, orbit.u0[2], 0.0, orbit.u0[3], 0.0]

"""
    get_orbit_length(orbit::TypeAPeriodicOrbit, type::Type{AbstractIndependantVariable})

Get the length of the periodic `orbit` in the desired variable.

# Arguments
- `orbit::TypeAPeriodicOrbit`: The periodic orbit.
- `type::Type{AbstractIndependantVariable}`: The type of length to return.

# Returns
- `length::Float64`: The length of the periodic orbit in the desired variable.
"""
get_orbit_length(orbit::AbstractPeriodicOrbit, ::Type{Time}) = orbit.P
get_orbit_length(orbit::AbstractPeriodicOrbit, ::Type{ArcLength}) = orbit.S

"""
    get_jacobi_integral(orbit::AbstractPeriodicOrbit)

Get the Jacobi integral of the periodic `orbit`.

# Arguments
- `orbit::AbstractPeriodicOrbit`: The periodic orbit.

# Returns
- `C::Float64`: The Jacobi integral of the periodic orbit.
"""
get_jacobi_integral(orbit::AbstractPeriodicOrbit) = jacobi_integral(get_full_initial_state(orbit), orbit.mu)

"""
    propagate_from_initial_conditions(
        orbit::AbstractPeriodicOrbit, τ1::AbstractFloat, PT::ParameterizationType,
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    )

Propagate the periodic `orbit` from the initial conditions to another state on the
`orbit`, as parameterized by `τ1` ∈ [0, 1] in terms of `PT::ParameterizationType`.

# Arguments
- `orbit::AbstractPeriodicOrbit`: The periodic orbit.
- `τ1::AbstractFloat`: The parameterization of the final state on the orbit.
- `PT::ParameterizationType`: The type of parameterization to use.

# Keyword Arguments
- `solver::OrdinaryDiffEq.AbstractODESolver`: The ODE solver to use for the integration.
- `reltol::Real`: The relative tolerance for the ODE solver.
- `abstol::Real`: The absolute tolerance for the ODE solver.

# Returns
- `x1::SVector{6,Float64}`: The final state on the periodic orbit.
"""
function propagate_from_initial_conditions(
    orbit::AbstractPeriodicOrbit, τ1::AbstractFloat, PT::ParameterizationType;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Check that τ ∈ [0, 1]
    if τ1 < 0.0 || τ1 > 1.0
        error("τ1 must be in [0, 1]")
    end

    # Form initial state and time span
    x0 = get_full_initial_state(orbit)
    tspan = (0.0, τ1*get_orbit_length(orbit, PT))

    # Propagate and return
    if τ1 == 0.0
        return x0
    elseif τ1 == 1.0
        return x0
    else
        return propagate_return_final_state(
            x0, tspan, orbit.mu, PT;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        )
    end
end

"""
    get_state_and_monodromy_matrix(orbit::AbstractPeriodicOrbit, τ1::AbstractFloat, PT::ParameterizationType)

Get the state and monodromy matrix of the periodic `orbit` at a give state
parameterized by `τ1` ∈ [0, 1] in terms of `PT`.

# Arguments
- `orbit::AbstractPeriodicOrbit`: The periodic orbit.
- `τ1::AbstractFloat`: The parameterization of the final state on the orbit.
- `PT::ParameterizationType`: The type of parameterization to use.
"""
function get_state_and_monodromy_matrix(
    orbit::AbstractPeriodicOrbit, τ1::AbstractFloat, PT::ParameterizationType;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    x0 = propagate_from_initial_conditions(
        orbit, τ1, PT;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
    M = propagate_return_final_stm(
        x0, (0.0, orbit.P), orbit.mu;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
    return x0, M
end

"""
    get_full_orbit(orbit::AbstractPeriodicOrbit)

Get the full orbit of the periodic `orbit`.

# Arguments
- `orbit::AbstractPeriodicOrbit`: The periodic orbit.
"""
function get_full_orbit(
    orbit::AbstractPeriodicOrbit;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    return propagate_return_all_states(
        get_full_initial_state(orbit),
        (0.0, orbit.P),
        orbit.mu,
        Time;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

function Base.show(
    io::IO, orbit::TypeAPeriodicOrbit,
)
    compact = get(io, :compact, false)
    if compact
        println(io, "TypeAPeriodicOrbit")
    else
        println(io, "TypeAPeriodicOrbit")
        println(io, "  Orbital period:  $(orbit.P)")
        println(io, "  Arc-length:      $(orbit.S)")
        println(io, "  N_cross:         $(orbit.N_cross)")
    end
end
