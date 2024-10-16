
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

    # The jacobi integral of the periodic orbit
    jacobi::Float64

    # The stability of the orbit
    stability::Float64

    # The period of the periodic orbit
    P::Float64

    # The arc length along the periodic orbit
    S::Float64

    # The number of crossings through x-z plane in half period
    N_cross::Int

    # The gravitational parameter
    mu::Float64

    # Units
    TU::Float64
    DU::Float64
end

"""
    TypeAPeriodicOrbit(
        rx, rz, vy, mu[, TU, DU];
        N_cross = 1,
        constraint = (:x_start_coordinate, 0.8),
        ode_solver = Vern9(),
        ode_reltol = 1e-12,
        ode_abstol = 1e-12,
        nl_solver  = SimpleTrustRegion(autodiff = nothing, nlsolve_update_rule = Val(true)),
    )
"""
function TypeAPeriodicOrbit(
    rx, rz, vy, mu, TU = NaN, DU = NaN;
    N_cross = 1,
    constraint = (:x_start_coordinate, rx),
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

    # Compute jacobi integral
    jacobi = jacobi_integral(SA[rx, 0.0, rz, 0.0, vy, 0.0], mu)

    # Compute the stability index
    M = propagate_return_final_stm(
        SA[rx, 0.0, rz, 0.0, vy, 0.0], (0.0, P), mu;
        solver = ode_solver,
        reltol = ode_reltol,
        abstol = ode_abstol,
    )
    vals, _ = eigen(M)
    idx_min = argmin(abs.(real(vals)))
    idx_max = argmax(abs.(real(vals)))
    stability = 0.5*(abs(vals[idx_min]) + abs(vals[idx_max]))

    return TypeAPeriodicOrbit(
        SA[rx,rz,vy], jacobi, stability, P, S, N_cross, mu, TU, DU,
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
        orbit.u0[1], orbit.u0[2], orbit.u0[3], orbit.mu, orbit.TU, orbit.DU;
        N_cross = orbit.N_cross,
        constraint = constraint,
        ode_solver = ode_solver,
        ode_reltol = ode_reltol,
        ode_abstol = ode_abstol,
        nl_solver  = nl_solver,
    )
end

"""
    GeneralPeriodicOrbit

A general representation of a periodic orbit.

# Fields
- `x0::SVector{6, Float64}`: The initial conditions for the periodic orbit.
- `P::Float64`: The period of the periodic orbit.
- `S::Float64`: The arc length along the periodic orbit.
- `mu::Float64`: The gravitational parameter of the CRTBP.
"""
struct GeneralPeriodicOrbit <: AbstractPeriodicOrbit
    # The initial conditions for the periodic orbit
    x0::SVector{6, Float64}

    # The jacobi integral of the periodic orbit
    jacobi::Float64

    # The stability of the orbit
    stability::Float64

    # The period of the periodic orbit
    P::Float64

    # The arc length along the periodic orbit
    S::Float64

    # The gravitational parameter
    mu::Float64

    # Units
    TU::Float64
    DU::Float64
end

function GeneralPeriodicOrbit(
    x0, P, mu;
    jacobi     = NaN,
    stability  = NaN,
    TU         = NaN,
    DU         = NaN,
    ode_solver = Vern9(),
    ode_reltol = 1e-14,
    ode_abstol = 1e-14,
)
    # Compute orbit arc length
    S = propagate_return_arc_length(
        x0, (0.0, P), mu;
        solver = ode_solver,
        reltol = ode_reltol,
        abstol = ode_abstol,
    )

    return GeneralPeriodicOrbit(x0, jacobi, stability, P, S, mu, TU, DU)
end

"""
    get_full_initial_state(orbit::AbstractPeriodicOrbit)

Get the full initial reference state of the periodic `orbit`.

# Arguments
- `orbit::AbstractPeriodicOrbit`: The periodic orbit.

# Returns
- `x0::SVector{6,Float64}`: The full initial state of the periodic orbit.
"""
get_full_initial_state(orbit::TypeAPeriodicOrbit) = SA[orbit.u0[1], 0.0, orbit.u0[2], 0.0, orbit.u0[3], 0.0]
get_full_initial_state(orbit::GeneralPeriodicOrbit) = orbit.x0

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

"""
    generate_periodic_orbit_cheb_interpolation(
        orbit::AbstractPeriodicOrbit, τ1_order::Int, PT::ParameterizationType;
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) -> FastChebInterpolation{FastChebInterp.ChebPoly{1, SVector{6, Float64}, Float64}}

Generate a Chebyshev interpolation of the periodic `orbit` in terms of `PT` employing
a Chebyshev polynomial of order `τ1_order`.

# Arguments
- `orbit::AbstractPeriodicOrbit`: The periodic orbit.
- `τ1_order::Int`: The order of the Chebyshev polynomial.
- `PT::ParameterizationType`: The type of parameterization to use (`Time` or `ArcLength`).

# Keyword Arguments
- `solver::OrdinaryDiffEq.AbstractODESolver`: The ODE solver to use for the integration.
- `reltol::Real`: The relative tolerance for the ODE solver.
- `abstol::Real`: The absolute tolerance for the ODE solver.

# Returns
- `FastChebInterpolation{FastChebInterp.ChebPoly{1, SVector{6, Float64}, Float64}}`: The Chebyshev interpolation.
"""
function generate_periodic_orbit_cheb_interpolation(
    orbit::AbstractPeriodicOrbit, τ1_order::Int, PT::ParameterizationType;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Get chebyshev nodes
    τ1s = chebpoints(τ1_order, 0.0, 1.0)

    # Compute state for each node
    states = Vector{SVector{6,Float64}}(undef, τ1_order + 1)
    for (i, τ1) in enumerate(τ1s)
        states[i] = propagate_from_initial_conditions(
            orbit, τ1, PT;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        )
    end

    # Return interpolation
    return FastChebInterpolation(chebinterp(states, 0.0, 1.0))
end

"""
    generate_periodic_orbit_cheb_approximation(
        orbit::AbstractPeriodicOrbit, τ1_npoints::Int, τ1_order::Int, PT::ParameterizationType;
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) -> FastChebInterpolation{FastChebInterp.ChebPoly{1, SVector{6, Float64}, Float64}}

Generate a Chebyshev approximation of the periodic `orbit` in terms of `PT` employing
`τ1_npoints` fit to a Chebyshev polynomial of order `τ1_order`.

# Arguments
- `orbit::AbstractPeriodicOrbit`: The periodic orbit.
- `τ1_npoints::Int`: The number of points to fit the Chebyshev polynomial.
- `τ1_order::Int`: The order of the Chebyshev polynomial.
- `PT::ParameterizationType`: The type of parameterization to use (`Time` or `ArcLength`).

# Keyword Arguments
- `solver::OrdinaryDiffEq.AbstractODESolver`: The ODE solver to use for the integration.
- `reltol::Real`: The relative tolerance for the ODE solver.
- `abstol::Real`: The absolute tolerance for the ODE solver.

# Returns
- `FastChebInterpolation{FastChebInterp.ChebPoly{1, SVector{6, Float64}, Float64}}`: The Chebyshev approximation.
"""
function generate_periodic_orbit_cheb_approximation(
    orbit::AbstractPeriodicOrbit,
    τ1_npoints::Int,
    τ1_order::Int,
    PT::ParameterizationType;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Get chebyshev nodes
    τ1s = chebpoints(τ1_npoints - 1, 0.0, 1.0)

    # Compute state for each node
    states = Vector{SVector{6,Float64}}(undef, τ1_npoints)
    for (i, τ1) in enumerate(τ1s)
        states[i] = propagate_from_initial_conditions(
            orbit, τ1, PT;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        )
    end

    # Return interpolation
    return FastChebInterpolation(chebregression(τ1s, states, 0.0, 1.0, τ1_order))
end

function Base.show(
    io::IO, ::MIME"text/plain", orbit::TypeAPeriodicOrbit,
)
    compact = get(io, :compact, false)
    if compact
        println(io, "TypeAPeriodicOrbit")
    else
        println(io, "TypeAPeriodicOrbit")
        if isnan(orbit.TU)
            println(io, "  Orbital period:  $(orbit.P) TU")
        else
            println(io, "  Orbital period:  $(orbit.P*orbit.TU/86400.0) days")
        end
        if isnan(orbit.DU)
            println(io, "  Arc-length:      $(orbit.S) DU")
        else
            println(io, "  Arc-length:      $(orbit.S*orbit.DU) km")
        end
        if !isnan(orbit.jacobi)
            println(io, "  Jacobi integral: $(orbit.jacobi)")
        end
        if !isnan(orbit.stability)
            println(io, "  Stability:       $(orbit.stability)")
        end
        println(io, "  N_cross:         $(orbit.N_cross)")
    end
end
function Base.show(
    io::IO, ::MIME"text/plain", orbit::GeneralPeriodicOrbit,
)
    compact = get(io, :compact, false)
    if compact
        println(io, "GeneralPeriodicOrbit")
    else
        println(io, "GeneralPeriodicOrbit")
        if isnan(orbit.TU)
            println(io, "  Orbital period:  $(orbit.P) TU")
        else
            println(io, "  Orbital period:  $(orbit.P*orbit.TU/86400.0) days")
        end
        if isnan(orbit.DU)
            println(io, "  Arc-length:      $(orbit.S) DU")
        else
            println(io, "  Arc-length:      $(orbit.S*orbit.DU) km")
        end
        if !isnan(orbit.jacobi)
            println(io, "  Jacobi integral: $(orbit.jacobi)")
        end
        if !isnan(orbit.stability)
            println(io, "  Stability:       $(orbit.stability)")
        end
    end
end
