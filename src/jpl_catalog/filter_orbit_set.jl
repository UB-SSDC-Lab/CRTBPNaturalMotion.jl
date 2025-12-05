
function get_filter_event_condition(orbit_set::OrbitSet, min_r1, max_r1, min_r2, max_r2)
    return (out, u, t, integrator) -> begin
        xpmu = u[1] + orbit_set.system.mass_ratio
        xpmum1 = xpmu - 1.0
        toDU = 1.0 / orbit_set.system.DU
        r1 = sqrt(xpmu*xpmu + u[2]*u[2] + u[3]*u[3])
        r2 = sqrt(xpmum1*xpmum1 + u[2]*u[2] + u[3]*u[3])
        out[1] = ifelse(isinf(min_r1), 1.0, r1 - min_r1*toDU)
        out[2] = ifelse(isinf(max_r1), 1.0, max_r1*toDU - r1)
        out[3] = ifelse(isinf(min_r2), 1.0, r2 - min_r2*toDU)
        out[4] = ifelse(isinf(max_r2), 1.0, max_r2*toDU - r2)
        return nothing
    end
end

"""
    filter_orbit_set(
        orbit_set::OrbitSet;
        min_r1      = -Inf,
        max_r1      = Inf,
        min_r2      = orbit_set.system.radius_secondary,
        max_r2      = Inf,
        ode_solver  = Vern9(),
        ode_reltol  = 1e-14,
        ode_abstol  = 1e-14,
    )

Filters a given `OrbitSet` with additional criteria that are more expensive to evaluate
than what we'd like to consider in `get_jpl_orbits()`

# Arguments
- `orbit_set::OrbitSet`: The orbit set you'd like to filter

# Keyword Arguments
- `min_r1::AbstractFloat = nothing`: The minimum allowable radius (in km) about the primary body.
- `max_r1::AbstractFloat = nothing`: The maximum allowable radius (in km) about the primary body.
- `min_r2::AbstractFloat = orbit_set.system.radius_secondary`: The minimum allowable radius (in km) about the secondary body.
- `max_r2::AbstractFloat = nothing`: The maximum allowable radius (in km) about the secondary body.
- `ode_solver = Vern9()`: The ODE solver to use for checking criteria along a full orbit.
- `ode_reltol = 1e-14`: Relative tolerance for the ODE solver. Default is `1e-14`.
- `ode_abstol = 1e-14`: Absolute tolerance for the ODE solver. Default is `1e-14`.
"""
function filter_orbit_set(
    orbit_set::OrbitSet;
    min_r1=(-Inf),
    max_r1=Inf,
    min_r2=orbit_set.system.radius_secondary,
    max_r2=Inf,
    ode_solver=Vern9(),
    ode_reltol=1e-14,
    ode_abstol=1e-14,
)
    # Get the ODE event condition for filtering
    filter_cond = get_filter_event_condition(orbit_set, min_r1, max_r1, min_r2, max_r2)

    # Construct callback
    filter_cb = VectorContinuousCallback(filter_cond, (integ, idx)->terminate!(integ), 4)

    # Loop through orbits and keep those that satisfy the criteria
    out = MVector{4,Float64}(undef)
    new_orbits = Vector{GeneralPeriodicOrbit}(undef, 0)
    for orbit in orbit_set
        # First verify if initial state satisfies criteria
        filter_cond(out, orbit.x0, 0.0, nothing)
        if any(SVector(out) .<= 0.0)
            continue
        end

        # Attempt to propagate for one full orbit
        ode_sol = solve(
            ODEProblem{false,SciMLBase.FullSpecialize}(
                natural_crtbp_eom, orbit.x0, (0.0, orbit.P), (orbit_set.system.mass_ratio,)
            ),
            ode_solver;
            callback=filter_cb,
            abstol=ode_abstol,
            reltol=ode_reltol,
            save_start=false,
            save_everystep=false,
        )

        # Check if we successfully propagated along the full orbit without terminating
        if isapprox(ode_sol.t[end], orbit.P)
            push!(new_orbits, orbit)
        end
    end

    return OrbitSet(new_orbits, orbit_set.system)
end
