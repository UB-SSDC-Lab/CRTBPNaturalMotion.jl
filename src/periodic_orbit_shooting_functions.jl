
"""
    typeA_shooting_function(u, p; constraint = (:x_start_coordinate, 0.0))

Shooting function for computing Type A periodic orbits in the CRTBP using the methodology
described by Restrepo and Russell in *A database of planar axisymmetric periodic orbits for
the solar system* [doi: 10.1007/s10569-018-9844-6](https://doi.org/10.1007/s10569-018-9844-6).

# Arguments
- `u::SVector{4,T}`: Vector of decision variables, [rx, rz, vy, P / 2], where P is half of
    the orbital period.
- `p::Tuple{Float64,}`: Tuple of parameters, (mu,), where mu is the mass parameter.

# Keyword Arguments
- `constraint::Tuple{Symbol,Real}`: Tuple of the constraint type and value corresponding
    to the final constraint considered to provide a fully defined system of nonlinear
    eqautions. The constraint type can be either `:x_start_coordinate` or `:jacobi_integral`.
    The constraint value is the value that the constraint should be equal to for the desired
    periodic orbit. Note that new constraints can be added with relative ease through
    modification of this function and its corresponding jacobian function
    `typeA_shooting_jacobian`.
- `solver::OrdinaryDiffEq.AbstractODESolver`: The solver to use for the integration.
- `reltol`: The relative tolerance for the solver.
- `abstol`: The absolute tolerance for the solver.

# Returns
- `SVector{4,T}`: Vector of residuals for the nonlinear system of equations.
"""
function typeA_shooting_function(
    u::SVector{4,T}, p;
    constraint  = (:x_start_coordinate, 0.0),
    solver      = Vern9(),
    reltol      = 1e-14,
    abstol      = 1e-14,
) where T
    # Grab the parameter
    mu = p[1]

    # Define spacecraft initial state
    x0 = SA[u[1], 0.0, u[2], 0.0, u[3], 0.0]

    # Grab the half-period
    half_period = u[4]

    # Propagate dynamics for half-period
    xf = propagate_return_final_state(
        x0, (0.0, half_period), mu;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )

    # Compute the extra constraint
    F4 = zero(eltype(u))
    if constraint[1] == :x_start_coordinate
        F4 += u[1] - constraint[2]
    elseif constraint[1] == :jacobi_integral
        F4 += jacobi_integral(x0, mu) - constraint[2]
    else
        error("Constraint type $(string(constraint[1])) not recognized.")
    end

    # Full constraint vector and return
    return SA[xf[2], xf[4], xf[6], F4]
end
function typeA_shooting_function!(
    du::AbstractVector{T}, u::AbstractVector{T}, p;
    constraint  = (:x_start_coordinate, 0.0),
    solver      = Vern9(),
    reltol      = 1e-14,
    abstol      = 1e-14,
) where T
    du .= typeA_shooting_function(
        SA[u[1], u[2], u[3], u[4]], p;
        constraint = constraint,
        solver     = solver,
        reltol     = reltol,
        abstol     = abstol,
    )
end


"""
    typeA_shooting_jacobian(u, p; constraint = (:x_start_coordinate, 0.0))

The Jacobian of the shooting function for computing Type A periodic orbits in the CRTBP
using the methodology described by Restrepo and Russell in *A database of planar
axisymmetric periodic orbits for the solar system* [doi: 10.1007/s10569-018-9844-6](https://doi.org/10.1007/s10569-018-9844-6).

# Arguments
- `u::SVector{4,T}`: Vector of decision variables, [rx, rz, vy, P / 2], where P is half of
    the orbital period.
- `p::Tuple{Float64,}`: Tuple of parameters, (mu,), where mu is the mass parameter.

# Keyword Arguments
- `constraint::Tuple{Symbol,Real}`: Tuple of the constraint type and value corresponding
    to the final constraint considered to provide a fully defined system of nonlinear
    eqautions. The constraint type can be either `:x_start_coordinate` or `:jacobi_integral`.
- `solver::OrdinaryDiffEq.AbstractODESolver`: The solver to use for the integration.
- `reltol`: The relative tolerance for the solver.
- `abstol`: The absolute tolerance for the solver.

# Returns
- `SMatrix{4,4,T,16}`: Jacobian of the nonlinear system of equations.
"""
function typeA_shooting_jacobian(
    u::SVector{4,T}, p;
    constraint  = (:x_start_coordinate, 0.0),
    solver      = Vern9(),
    reltol      = 1e-14,
    abstol      = 1e-14,
) where T
    # Grab the parameter
    mu = p[1]

    # Define spacecraft initial state
    x0 = SA[u[1], 0.0, u[2], 0.0, u[3], 0.0]

    # Grab the half-period
    half_period = u[4]

    # Propagate dynamics for half-period to get STM
    xf, STM = propagate_return_final_state_and_stm(
        x0, (0.0, half_period), mu;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )

    # Get the dynamics at the final state
    dxf = natural_crtbp_eom(xf, (mu,), half_period)

    # Compute the extra constraint jacobian
    if constraint[1] == :x_start_coordinate
        J4 = SA[1.0, 0.0, 0.0]
    elseif constraint[1] == :jacobi_integral
        dC_dx0 = jacobi_integral_gradient(x0, mu)
        J4 = dC_dx0[SA[1,3,5]]
    else
        error("Constraint type $(string(constraint[1])) not recognized.")
    end

    # Compute and return Jacobian
    return SA[
        STM[2,1] STM[2,3] STM[2,5] dxf[2];
        STM[4,1] STM[4,3] STM[4,5] dxf[4];
        STM[6,1] STM[6,3] STM[6,5] dxf[6];
           J4[1]    J4[2]    J4[3]    0.0;
    ]
end
function typeA_shooting_jacobian!(
    J::AbstractMatrix{T}, u::AbstractVector{T}, p;
    constraint  = (:x_start_coordinate, 0.0),
    solver      = Vern9(),
    reltol      = 1e-14,
    abstol      = 1e-14,
) where T
    J .= typeA_shooting_jacobian(
        SA[u[1], u[2], u[3], u[4]], p;
        constraint = constraint,
        solver     = solver,
        reltol     = reltol,
        abstol     = abstol,
    )
end
