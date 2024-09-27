@inline synodic_frame_xz_plane_crossing_condition(u::SVector,t,integ) = u[2]
@inline synodic_frame_xz_plane_crossing_condition(u::SMatrix,t,integ) = u[2,1]

"""
    correct_typeA_initial_conditions(
        rx_guess, rz_guess, vy_guess, mu;
        N_cross     = 1,
        constraint  = (:x_start_coordinate, 0.0),
        ode_solver  = Vern9(),
        ode_reltol  = 1e-14,
        ode_abstol  = 1e-14,
        nl_solver   = SimpleTrustRegion(autodiff = nothing, nlsolve_update_rule = Val(true)),
    )

Correct the initial conditions of a type-A periodic orbit using a forward shooting method.

# Arguments
- `rx_guess::Real`: Initial guess for the x-coordinate of the initial state.
- `rz_guess::Real`: Initial guess for the z-coordinate of the initial state.
- `vy_guess::Real`: Initial guess for the y-velocity of the initial state.
- `mu::Real`: Gravitational parameter of the CRTBP.

# Keyword Arguments
- `N_cross::Int = 1`: Number of crossings of the xz-plane to consider in the half period.
- `constraint::Tuple{Symbol,Real}`: Tuple of the constraint type and value corresponding
    to the final constraint considered to provide a fully defined system of nonlinear
    eqautions. The constraint type can be either `:x_start_coordinate` or `:jacobi_integral`.
    The constraint value is the value that the constraint should be equal to for the desired
    periodic orbit.
- `ode_solver::OrdinaryDiffEq.AbstractODESolver`: The ODE solver to use for the integration.
- `ode_reltol::Real`: The relative tolerance for the ODE solver.
- `ode_abstol::Real`: The absolute tolerance for the ODE solver.
- `nl_solver::NonlinearSolve.AbstractNLsolveSolver`: The nonlinear solver to use for the
    root-finding problem. Note, there seems to be a strange bug with the NonlinearSolve.jl
    `TrustRegion` and `NewtonRaphson` solvers when using out-of-place Jacobian defined
    to return an SMatrix. Therefore, it is recommened to use either the default
    `SimpleTrustRegion` solver or the `SimpleNewtonRaphson` solver (which are specifically
    tailored towards StaticArrays.jl arrays).
"""
function correct_typeA_initial_conditions(
    rx_guess, rz_guess, vy_guess, mu;
    N_cross     = 1,
    constraint  = (:x_start_coordinate, rx_guess),
    ode_solver  = Vern9(),
    ode_reltol  = 1e-14,
    ode_abstol  = 1e-14,
    nl_solver   = TrustRegion(autodiff = nothing),
)
    # Define initial state from guess
    x0_guess = SA[rx_guess, 0.0, rz_guess, 0.0, vy_guess, 0.0]

    # Determine guess for half period
    xi = x0_guess
    half_P = 0.0
    for _ in 1:N_cross
        ode_sol = propagate_return_all_states(
            xi, (0.0, Inf), mu, synodic_frame_xz_plane_crossing_condition;
            solver = ode_solver,
            reltol = ode_reltol,
            abstol = ode_abstol,
        )
        half_P += ode_sol.t[end]
        xi = ode_sol.u[end]
    end

    # Construct shooting function and jacobian
    fun = let constraint = constraint
        (du,u,p) -> typeA_shooting_function!(du, u, p, constraint = constraint)
    end
    jac = let constraint = constraint
        (J,u,p) -> typeA_shooting_jacobian!(J, u, p, constraint = constraint)
    end
    nlfun = NonlinearFunction{true, SciMLBase.FullSpecialize}(fun; jac = jac)
    nlp = NonlinearProblem{true}(nlfun, [rx_guess,rz_guess,vy_guess,half_P], (mu,))

    # Solve and return solution with flag indicating success/failure (true/false)
    nlsol = solve(nlp, nl_solver)
    return (
        SA[nlsol.u[1], nlsol.u[2], nlsol.u[3], 2.0*nlsol.u[4]],
        SciMLBase.successful_retcode(nlsol)
    )
end
