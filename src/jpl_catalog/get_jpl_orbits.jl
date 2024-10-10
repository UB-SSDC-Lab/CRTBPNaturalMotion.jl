# ===== Struct for returning data to user
struct System
    name::String
    mass_ratio::Float64
    radius_secondary::Float64
    L1::SVector{3, Float64}
    L2::SVector{3, Float64}
    L3::SVector{3, Float64}
    L4::SVector{3, Float64}
    L5::SVector{3, Float64}
    LU::Float64
    TU::Float64

    function System(json_sys::JSONSystem)
        new(
            json_sys.name,
            json_sys.mass_ratio,
            json_sys.radius_secondary,
            SA[json_sys.L1[1], json_sys.L1[2], json_sys.L1[3]],
            SA[json_sys.L2[1], json_sys.L2[2], json_sys.L2[3]],
            SA[json_sys.L3[1], json_sys.L3[2], json_sys.L3[3]],
            SA[json_sys.L4[1], json_sys.L4[2], json_sys.L4[3]],
            SA[json_sys.L5[1], json_sys.L5[2], json_sys.L5[3]],
            json_sys.lunit,
            json_sys.tunit,
        )
    end
end

struct Limits
    jacobi::SVector{2, Float64}
    period::SVector{2, Float64}
    stability::SVector{2, Float64}

    function Limits(json_limits::JSONLimits)
        new(
            SA[json_limits.jacobi[1], json_limits.jacobi[2]],
            SA[json_limits.period[1], json_limits.period[2]],
            SA[json_limits.stability[1], json_limits.stability[2]],
        )
    end
    function Limits(
        jacobi::SVector{2, Float64},
        period::SVector{2, Float64},
        stability::SVector{2, Float64},
    )
        new(jacobi, period, stability)
    end
end

struct OrbitSet <: AbstractVector{GeneralPeriodicOrbit}
    # Vector of orbits
    orbits::Vector{GeneralPeriodicOrbit}

    # System data
    system::System

    # Limits data
    limits::Limits

    function OrbitSet(
        orbits::Vector{GeneralPeriodicOrbit},
        system::JSONSystem,
    )
        # Seems to be some issues with the returned limits, so we'll just compute them instead
        minperiod, maxperiod = Inf, -Inf
        minjacobi, maxjacobi = Inf, -Inf
        minstab, maxstab     = Inf, -Inf
        for orbit in orbits
            if orbit.P < minperiod
                minperiod = orbit.P
            end
            if orbit.P > maxperiod
                maxperiod = orbit.P
            end
            if orbit.jacobi < minjacobi
                minjacobi = orbit.jacobi
            end
            if orbit.jacobi > maxjacobi
                maxjacobi = orbit.jacobi
            end
            if orbit.stability < minstab
                minstab = orbit.stability
            end
            if orbit.stability > maxstab
                maxstab = orbit.stability
            end
        end
        new(
            orbits, System(system),
            Limits(
                SA[minjacobi, maxjacobi],
                SA[minperiod, maxperiod],
                SA[minstab, maxstab],
            ),
        )
    end
    function OrbitSet(
        orbits::Vector{GeneralPeriodicOrbit},
        system::System,
    )
        # Seems to be some issues with the returned limits, so we'll just compute them instead
        minperiod, maxperiod = Inf, -Inf
        minjacobi, maxjacobi = Inf, -Inf
        minstab, maxstab     = Inf, -Inf
        for orbit in orbits
            if orbit.P < minperiod
                minperiod = orbit.P
            end
            if orbit.P > maxperiod
                maxperiod = orbit.P
            end
            if orbit.jacobi < minjacobi
                minjacobi = orbit.jacobi
            end
            if orbit.jacobi > maxjacobi
                maxjacobi = orbit.jacobi
            end
            if orbit.stability < minstab
                minstab = orbit.stability
            end
            if orbit.stability > maxstab
                maxstab = orbit.stability
            end
        end
        new(
            orbits, system,
            Limits(
                SA[minjacobi, maxjacobi],
                SA[minperiod, maxperiod],
                SA[minstab, maxstab],
            ),
        )
    end
end

# Define interface
Base.size(os::OrbitSet) = size(os.orbits)
Base.getindex(os::OrbitSet, i::Int) = os.orbits[i]
Base.getindex(os::OrbitSet, rng::AbstractVector{Int}) = OrbitSet(os.orbits[rng], os.system)
Base.IndexStyle(::Type{OrbitSet}) = Base.IndexStyle(Vector{GeneralPeriodicOrbit})
Base.setindex!(os::OrbitSet, o::GeneralPeriodicOrbit, i::Int) = setindex!(os.orbits, o, i)
Base.length(os::OrbitSet) = length(os.orbits)

# Utility functions
function minimum_period_orbit(os::OrbitSet)
    min_period = Inf
    idx = 0
    for (i,orbit) in enumerate(os.orbits)
        if orbit.P < min_period
            min_period = orbit.P
            idx = i
        end
    end
    return os[idx]
end
function maximum_period_orbit(os::OrbitSet)
    max_period = -Inf
    idx = 0
    for (i,orbit) in enumerate(os.orbits)
        if orbit.P > max_period
            max_period = orbit.P
            idx = i
        end
    end
    return os[idx]
end

# Define pretty printing
function Base.show(io::IO, ::MIME"text/plain", os::OrbitSet)
    compact = get(io, :compact, false)
    if compact
        println(io, "OrbitSet with $(length(os)) orbits")
    else
        period_min_days = os.limits.period[1] * os.system.TU / 86400.0
        period_max_days = os.limits.period[2] * os.system.TU / 86400.0
        println(io, "OrbitSet with $(length(os)) orbits")
        printfmt(io, "\tSystem:            {1:s}\n", os.system.name)
        printfmt(io, "\tPeriod range:    {1:8.4f} to {2:.4f} days\n", period_min_days, period_max_days)
        printfmt(io, "\tJacobi range:    {1:8.4f} to {2:.4f}\n", os.limits.jacobi[1], os.limits.jacobi[2])
        printfmt(io, "\tStability range: {1:8.4f} to {2:.4f}\n", os.limits.stability[1], os.limits.stability[2])
    end
end

"""
    get_jpl_orbits(; sys = "earth-moon", family = "halo", libr = 1, branch = "N",
                     periodmin = nothing, periodmax = nothing, periodunits = nothing,
                     jacobimin = nothing, jacobimax = nothing, stabmin = nothing, stabmax = nothing,
                     ode_solver = Vern9(), ode_reltol = 1e-14, ode_abstol = 1e-14)

Get periodic orbits from the JPL CRTBP Periodic Orbit API. See JPL's [web application](https://ssd.jpl.nasa.gov/tools/periodic_orbits.html)
and [API documentation](https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html) for details.

# Keyword Arguments
- `sys::String = "earth-moon"`: The system to get orbits for. Options include:
    - `"earth-moon"`
    - `"mars-phobos"`
    - `"sun-earth"`
    - `"sun-mars"`
    - `"jupiter-europa"`
    - `"saturn-enceladus"`
    - `"saturn-titan"`
- `family::String = "halo"`: The family of orbits to get. Options include:
    - `"halo"`
    - `"vertical"`
    - `"axial"`
    - `"lyapunov"`
    - `"longp"`
    - `"short"`
    - `"butterfly"`
    - `"dragonfly"`
    - `"resonant"`
    - `"dro"`
    - `"dpo"`
    - `"lpo"`
- `libr::Int = 1`: The libration point to get orbits around. Default is `1`. Requred for
    `"halo"`, `"lyapunov"` (`1`, `2`, or `3`), `"longp"`, `"short"` (`4` or `5`), and `"axial"`,
    `"vertical"` (`1`, `2`, `3`, `4`, or `5`).
- `branch::String = "N"`: The branch of the family to get. Required for `"halo"`,
    `"dragonfly"`, `"butterfly"` (`"N"` of `"S"`), `"lpo"` (`"E"` or `"W"`), and
    `"resonant"` (e.g., `"12"` for 1:2 resonant orbits).
- `periodmin::Union{Float64,Nothing} = nothing`: Minimum period of orbits to get. Default is `nothing`.
- `periodmax::Union{Float64,Nothing} = nothing`: Maximum period of orbits to get. Default is `nothing`.
- `periodunits::Union{String,Nothing} = nothing`: Units of period. Options include:
    - `"s"`: seconds
    - `"h"`: hours
    - `"d"`: days
    - `"TU"`: dimensionless time units
- `jacobimin::Union{Float64,Nothing} = nothing`: Minimum Jacobi constant of orbits to get. Default is `nothing`.
- `jacobimax::Union{Float64,Nothing} = nothing`: Maximum Jacobi constant of orbits to get. Default is `nothing`.
- `stabmin::Union{Float64,Nothing} = nothing`: Minimum stability of orbits to get. Default is `nothing`.
- `stabmax::Union{Float64,Nothing} = nothing`: Maximum stability of orbits to get. Default is `nothing`.
- `ode_solver = Vern9()`: The ODE solver to use for computing the orbits. Default is `Vern9()`. This is
    only employed for computing the arc-length along the orbit.
- `ode_reltol = 1e-14`: Relative tolerance for the ODE solver. Default is `1e-14`.
- `ode_abstol = 1e-14`: Absolute tolerance for the ODE solver. Default is `1e-14`.

# Returns
- `OrbitSet`: A set of periodic orbits that meet the query criteria.
"""
function get_jpl_orbits(;
    sys::String                        = "earth-moon",
    family::String                     = "halo",
    libr::Int                          = 1,
    branch::String                     = "N",
    periodmin::Union{Float64,Nothing}  = nothing,
    periodmax::Union{Float64,Nothing}  = nothing,
    periodunits::Union{String,Nothing} = nothing,
    jacobimin::Union{Float64,Nothing}  = nothing,
    jacobimax::Union{Float64,Nothing}  = nothing,
    stabmin::Union{Float64,Nothing}    = nothing,
    stabmax::Union{Float64,Nothing}    = nothing,
    ode_solver                         = Vern9(),
    ode_reltol                         = 1e-14,
    ode_abstol                         = 1e-14,
)
    return get_jpl_orbits(
        construct_query(;
            sys         = sys,
            family      = family,
            libr        = libr,
            branch      = branch,
            periodmin   = periodmin,
            periodmax   = periodmax,
            periodunits = periodunits,
            jacobimin   = jacobimin,
            jacobimax   = jacobimax,
            stabmin     = stabmin,
            stabmax     = stabmax,
        );
        ode_solver  = ode_solver,
        ode_reltol  = ode_reltol,
        ode_abstol  = ode_abstol,
    )
end

"""


Get periodic orbits from the JPL CRTBP Periodic Orbit API. See JPL's [web application](https://ssd.jpl.nasa.gov/tools/periodic_orbits.html)
and [API documentation](https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html) for details.

# Arguments
- `query::String`: The query string to send to the JPL API. See JPL's [API documentation](https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html)
    for details on manually constructing an HTTP Request query.

# Keyword Arguments
- `ode_solver = Vern9()`: The ODE solver to use for computing the orbits. Default is `Vern9()`. This is
    only employed for computing the arc-length along the orbit.
- `ode_reltol = 1e-14`: Relative tolerance for the ODE solver. Default is `1e-14`.
- `ode_abstol = 1e-14`: Absolute tolerance for the ODE solver. Default is `1e-14`.

# Returns
- `OrbitSet`: A set of periodic orbits that meet the query criteria.
"""
function get_jpl_orbits(
    query::String;
    ode_solver = Vern9(),
    ode_reltol = 1e-14,
    ode_abstol = 1e-14,
)
    # Request data with query
    orbit_data = request_data(query)

    # Get mass ratio and units
    mu = orbit_data.system.mass_ratio
    TU = orbit_data.system.tunit
    LU = orbit_data.system.lunit

    # Allocate memory for orbits and compute
    orbits = Vector{GeneralPeriodicOrbit}(undef, orbit_data.count)
    for i in 1:orbit_data.count
        # Form initial state
        x0 = SA[
            orbit_data.data[i][1], orbit_data.data[i][2], orbit_data.data[i][3],
            orbit_data.data[i][4], orbit_data.data[i][5], orbit_data.data[i][6],
        ]

        # Get period, jacobi, and stability
        P           = orbit_data.data[i][8]
        jacobi      = orbit_data.data[i][7]
        stability   = orbit_data.data[i][9]

        # Compute orbit
        orbits[i] = GeneralPeriodicOrbit(
            x0, P, mu;
            TU         = TU,
            LU         = LU,
            jacobi     = jacobi,
            stability  = stability,
            ode_solver = ode_solver,
            ode_reltol = ode_reltol,
            ode_abstol = ode_abstol,
        )
    end
    return OrbitSet(orbits, orbit_data.system)
end
