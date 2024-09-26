"""
    propagate_return_all_states(
        x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu [, iv::Type{AbstractIndependantVariable}];
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) where {TX,TT}

Propagate the state `x0` over the time span `tspan` and return the full ODE solution. The
final argument can be specified to set the desired independant variable (i.e., `Time` or
`ArcLength`).

# Arguments
- `x0::SVector{6,T}`: The initial state with x0 = [r0; v0].
- `tspan::Tuple{T,T}`: The span of the independant variable forwhich to solve the ODE.
- `mu::Real`: The mass parameter of the three-body system (`mu = m2 / (m1 + m2)`).
- `iv::Type{AbstractIndependantVariable}`: The desired independant variable (Time if
    unspecified).

# Keywords
- `solver::OrdinaryDiffEq.AbstractODESolver`: The solver to use for the integration.
- `reltol`: The relative tolerance for the solver.
- `abstol`: The absolute tolerance for the solver.

# Returns
- `ODESolution`: The full solution to the ODE problem.
"""
@inline function propagate_return_all_states(
    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu, iv::Type{Time};
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom,
            x0,
            tspan,
            (mu,),
        ),
        solver;
        abstol = abstol,
        reltol = reltol,
    )
end

@inline function propagate_return_all_states(
    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return propagate_return_all_states(
        x0, tspan, mu, Time;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

@inline function propagate_return_all_states(
    x0::SVector{6,TX}, sspan::Tuple{TT,TT}, mu, iv::Type{ArcLength};
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom_with_independant_arclen,
            x0,
            sspan,
            (mu,),
        ),
        solver;
        abstol = abstol,
        reltol = reltol,
    )
end

"""
    propagate_return_all_states(
        x0::SVector{6,TX},
        tsteps::Union{StepRangeLen,AbstractVector},
        mu,
        iv::Type{AbstractIndependantVariable};
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) where TX

Propagate the state `x0` from `tsteps[1]` to `tsteps[end]`, returning the solution at each
time in `tsteps`. The final argument can be specified to set the desired independant variable
(i.e., `Time` or `ArcLength`).

# Arguments
- `x0::SVector{6,TX}`: The initial state with x0 = [r0; v0].
- `tsteps::Union{StepRangeLen,AbstractVector}`: The time steps at which a state x(t) is desired.
- `mu::Real`: The mass parameter of the three-body system (`mu = m2 / (m1 + m2)`).
- `iv::Type{AbstractIndependantVariable}`: The desired independant variable.

# Keywords
- `solver::OrdinaryDiffEq.AbstractODESolver`: The solver to use for the integration.
- `reltol`: The relative tolerance for the solver.
- `abstol`: The absolute tolerance for the solver.

# Returns
- `Vector{SVector{6,TX}}`: A vector of states at each time in `tsteps`.
"""
@inline function propagate_return_all_states(
    x0::SVector{6,TX}, tsteps::Union{StepRangeLen,AbstractVector}, mu, iv::Type{Time};
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where TX
    ode_sol = solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom,
            x0,
            (tsteps[0], tsteps[end]),
            (mu,),
        ),
        solver;
        tstops = tsteps,
        abstol = abstol,
        reltol = reltol,
    )
    return ode_sol(tsteps)
end

@inline function propagate_return_all_states(
    x0::SVector{6,TX}, ssteps::Union{StepRangeLen,AbstractVector}, mu, ::Type{ArcLength};
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where TX
    ode_sol = solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom_with_independant_arclen,
            x0,
            (ssteps[1], ssteps[end]),
            (mu,),
        ),
        solver;
        tstops = ssteps,
        abstol = abstol,
        reltol = reltol,
    )
    return ode_sol(ssteps)
end

"""
    propagate_return_all_states(
        x0::SVector{6,TX},
        tspan::Tuple{TT,TT},
        mu,
        term_cond::Function[, iv::Type{AbstractIndependantVariable}];
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) where {TX,TT}

Propagate the state `x0` over the time span `tspan` and return the full ODE solution. At the
first point at which `term_cond` is zero (if any), the integration will terminate. The
final argument can be specified to set the desired independant variable (i.e., `Time` or
`ArcLength`).

# Arguments
- `x0::SVector{6,TX}`: The initial state with x0 = [r0; v0].
- `tspan::Tuple{TT,TT}`: The span of the independant variable forwhich to solve the ODE.
- `mu::Real`: The mass parameter of the three-body system (`mu = m2 / (m1 + m2)`).
- `term_cond::Function`: The DifferentialEquations.jl continuous callback condition function
    that will terminate the integration when it is zero. Function should be of the form
    `term_cond(x,t,integ)`, where `x` is the state, `t` is the independant variable, and
    `integ` is the integrator struct (see [DifferentialEquations.jl callback documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)).
- `iv::Type{AbstractIndependantVariable}`: The desired independant variable (Time if
    unspecified).

# Keywords
- `solver::OrdinaryDiffEq.AbstractODESolver`: The solver to use for the integration.
- `reltol`: The relative tolerance for the solver.
- `abstol`: The absolute tolerance for the solver.

# Returns
- `ODESolution`: The full solution to the ODE problem.
"""
@inline function propagate_return_all_states(
    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu, term_cond::Function, iv::Type{Time};
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom,
            x0,
            tspan,
            (mu,),
        ),
        solver;
        callback = ContinuousCallback(term_cond, integ -> terminate!(integ)),
        abstol = abstol,
        reltol = reltol,
    )
end

@inline function propagate_return_all_states(
    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu, term_cond::Function;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return propagate_return_all_states(
        x0, tspan, mu, term_cond, Time;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

@inline function propagate_return_all_states(
    x0::SVector{6,TX}, sspan::Tuple{TT,TT}, mu, term_cond::Function, iv::Type{ArcLength};
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom_with_independant_arclen,
            x0,
            sspan,
            (mu,),
        ),
        solver;
        callback = ContinuousCallback(term_cond, integ -> terminate!(integ)),
        abstol = abstol,
        reltol = reltol,
    )
end

"""
    propagate_return_all_states(
        x0::SVector{6,TX},
        tspan::Tuple{TT,TT},
        mu,
        cond::Function,
        affect::Function,
        iv::Type{AbstractIndependantVariable};
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) where {TX,TT}

Propagate the state `x0` over the time span `tspan` and return the full ODE solution. This
is done employing a DifferentialEquations.jl ContinuousCallback defined with `cond` and
`affect`. The final argument can be specified to set the desired independant variable
(i.e., `Time` or `ArcLength`).

# Arguments
- `x0::SVector{6,TX}`: The initial state with x0 = [r0; v0].
- `tspan::Tuple{TT,TT}`: The span of the independant variable forwhich to solve the ODE.
- `mu::Real`: The mass parameter of the three-body system (`mu = m2 / (m1 + m2)`).
- `cond::Function`: The DifferentialEquations.jl continuous callback condition function
    that will apply `affect` to the integrator struct (see
    [DifferentialEquations.jl callback documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)).
- `affect::Function`: The DifferentialEquations.jl continuous callback affect function.
- `iv::Type{AbstractIndependantVariable}`: The desired independant variable (Time if
    unspecified).

# Keywords
- `solver::OrdinaryDiffEq.AbstractODESolver`: The solver to use for the integration.
- `reltol`: The relative tolerance for the solver.
- `abstol`: The absolute tolerance for the solver.

# Returns
- `ODESolution`: The full solution to the ODE problem.
"""
@inline function propagate_return_all_states(
    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu,
    cond::Function, affect::Function,
    iv::Type{Time};
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom,
            x0,
            tspan,
            (mu,),
        ),
        solver;
        callback = ContinuousCallback(cond, affect),
        abstol = abstol,
        reltol = reltol,
    )
end

@inline function propagate_return_all_states(
    x0::SVector{6,TX}, sspan::Tuple{TT,TT}, mu,
    cond::Function, affect::Function,
    iv::Type{ArcLength};
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom_with_independant_arclen,
            x0,
            sspan,
            (mu,),
        ),
        solver;
        callback = ContinuousCallback(cond, affect),
        abstol = abstol,
        reltol = reltol,
    )
end

"""
    propagate_return_final_state(
        x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu [, iv::Type{AbstractIndependantVariable}];
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) where {TX,TT}

Propagate the state `x0` over the time span `tspan` and return the final state. The
final argument can be specified to set the desired independant variable (i.e., `Time` or
`ArcLength`).

# arguments
- `x0::SVector{6,TX}`: The initial state with x0 = [r0; v0].
- `tspan::Tuple{TT,TT}`: The span of the independant variable forwhich to solve the ode.
- `mu::Real`: The mass parameter of the three-body system (`mu = m2 / (m1 + m2)`).
- `iv::Type{AbstractIndependantVariable}`: The desired independant variable (Time if
    unspecified).

# keywords
- `solver::OrdinaryDiffEq.AbstractODESolver`: The solver to use for the integration.
- `reltol`: The relative tolerance for the solver.
- `abstol`: The absolute tolerance for the solver.

# returns
- `SVector{6,TX}`: The final state of the ode solution.
"""
@inline function propagate_return_final_state(
    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu, ::Type{Time};
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom,
            x0,
            tspan,
            (mu,),
        ),
        solver;
        abstol = abstol,
        reltol = reltol,
        save_start = false,
        save_everystep = false,
    )[end]
end

@inline function propagate_return_final_state(
    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return propagate_return_final_state(
        x0, tspan, mu, Time;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

@inline function propagate_return_final_state(
    x0::SVector{6,TX}, sspan::Tuple{TT,TT}, mu, ::Type{ArcLength};
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom_with_independant_arclen,
            x0,
            sspan,
            (mu,),
        ),
        solver;
        abstol = abstol,
        reltol = reltol,
        save_start = false,
        save_everystep = false,
    )[end]
end

"""
    propagate_return_final_state(
        x0::SVector{6,TX},
        tspan::Tuple{TT,TT},
        mu,
        term_cond::Function [, iv::Type{AbstractIndependantVariable}];
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) where {TX,TT}

Propagate the state `x0` over the time span `tspan` and return the final state. At the
first point at which `term_cond` is zero (if any), the integration will terminate. The
final argument can be specified to set the desired independant variable (i.e., `Time` or
`ArcLength`).

# Arguments
- `x0::SVector{6,TX}`: The initial state with x0 = [r0; v0].
- `tspan::Tuple{TT,TT}`: The span of the independant variable forwhich to solve the ODE.
- `mu::Real`: The mass parameter of the three-body system (`mu = m2 / (m1 + m2)`).
- `term_cond::Function`: The DifferentialEquations.jl continuous callback condition function
    that will terminate the integration when it is zero. Function should be of the form
    `term_cond(x,t,integ)`, where `x` is the state, `t` is the independant variable, and
    `integ` is the integrator struct (see [DifferentialEquations.jl callback documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/callback_functions/)).
- `iv::Type{AbstractIndependantVariable}`: The desired independant variable (Time if
    unspecified).

# Keywords
- `solver::OrdinaryDiffEq.AbstractODESolver`: The solver to use for the integration.
- `reltol`: The relative tolerance for the solver.
- `abstol`: The absolute tolerance for the solver.

# Returns
- `SVector{6,TX}`: The final state of the ODE solution.
"""
@inline function propagate_return_final_state(
    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu, term_cond::Function, ::Type{Time};
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom,
            x0,
            tspan,
            (mu,),
        ),
        solver;
        callback = ContinuousCallback(term_cond, integ -> terminate!(integ)),
        abstol = abstol,
        reltol = reltol,
        save_start = false,
        save_everystep = false,
    )[end]
end

@inline function propagate_return_final_state(
    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu, term_cond::Function;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return propagate_return_final_state(
        x0, tspan, mu, term_cond, Time;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

@inline function propagate_return_final_state(
    x0::SVector{6,TX}, sspan::Tuple{TT,TT}, mu, term_cond::Function, ::Type{ArcLength};
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    return solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom_with_independant_arclen,
            x0,
            sspan,
            (mu,),
        ),
        solver;
        callback = ContinuousCallback(term_cond, integ -> terminate!(integ)),
        abstol = abstol,
        reltol = reltol,
        save_start = false,
        save_everystep = false,
    )[end]
end

"""
    propagate_return_arc_length(
        x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu;
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) where {TX,TT}

Propagate the state `x0` over the time span `tspan` and return the arclength along the
trajectory. The final argument can be specified to set the desired independant variable
(i.e., `Time` or `ArcLength`).

# arguments
- `x0::SVector{6,TX}`: the initial state with x0 = [r0; v0].
- `tspan::Tuple{TT,TT}`: the time span forwhich to solve the ode.
- `mu::Real`: the mass parameter of the three-body system (`mu = m2 / (m1 + m2)`).

# keywords
- `solver::OrdinaryDiffEq.AbstractODESolver`: The solver to use for the integration.
- `reltol`: The relative tolerance for the solver.
- `abstol`: The absolute tolerance for the solver.

# returns
- `TX`: The arc-length along the curve
"""
@inline function propagate_return_arc_length(
    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    # Create state with initial arc-length
    x0s = SA[x0[1], x0[2], x0[3], x0[4], x0[5], x0[6], 0.0]

    # Propagate
    xfs = solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom_with_arclen,
            x0s,
            tspan,
            (mu,),
        ),
        solver;
        abstol = abstol,
        reltol = reltol,
        save_start = false,
        save_everystep = false,
    )[end]

    return xfs[7]
end

"""
    propagate_return_time(
        x0::SVector{6,TX}, sspan::Tuple{TT,TT}, mu;
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) where {TX,TT}

Propagate the state `x0` over the arc-length span `sspan` and return the arclength along the
trajectory. The final argument can be specified to set the desired independant variable
(i.e., `Time` or `ArcLength`).

# arguments
- `x0::SVector{6,TX}`: the initial state with x0 = [r0; v0].
- `sspan::Tuple{TT,TT}`: the arc-length span forwhich to solve the ode.
- `mu::Real`: the mass parameter of the three-body system (`mu = m2 / (m1 + m2)`).

# keywords
- `solver::OrdinaryDiffEq.AbstractODESolver`: The solver to use for the integration.
- `reltol`: The relative tolerance for the solver.
- `abstol`: The absolute tolerance for the solver.

# returns
- `TX`: The time along the curve
"""
@inline function propagate_return_time(
    x0::SVector{6,TX}, sspan::Tuple{TT,TT}, mu;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    # Create state with initial time
    x0t = SA[x0[1], x0[2], x0[3], x0[4], x0[5], x0[6], 0.0]

    xfs = solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_and_time_eom_with_independant_arclen,
            x0t,
            sspan,
            (mu,),
        ),
        solver;
        abstol = abstol,
        reltol = reltol,
        save_start = false,
        save_everystep = false,
    )[end]

    return xfs[7]
end

"""
    propagate_return_final_stm(
        x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu;
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) where {TX,TT}

Propagate the state `x0` over the time span `tspan` and return the final state transition
matrix (STM).

# arguments
- `x0::SVector{6,TX}`: The initial state with x0 = [r0; v0].
- `tspan::Tuple{TT,TT}`: The span of the independant variable forwhich to solve the ode.
- `mu::Real`: The mass parameter of the three-body system (`mu = m2 / (m1 + m2)`).

# keywords
- `solver::OrdinaryDiffEq.AbstractODESolver`: The solver to use for the integration.
- `reltol`: The relative tolerance for the solver.
- `abstol`: The absolute tolerance for the solver.

# returns
- `SMatrix{6,6,TX,36}`: The stm for the full trajectory
"""
@inline function propagate_return_final_stm(
    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    # Form initial state with STM
    z_start = initial_state_with_stm(x0)

    # Propagate
    z_stop = solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom_with_stm,
            z_start,
            tspan,
            (mu,),
        ),
        solver;
        abstol = abstol,
        reltol = reltol,
        save_start = false,
        save_everystep = false,
    )[end]

    return get_stm(z_stop)
end

"""
    propagate_return_final_state_and_stm(
        x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu;
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) where {TX,TT}

Propagate the state `x0` over the time span `tspan` and return the final state and the
state transition matrix (STM).

# arguments
- `x0::SVector{6,TX}`: The initial state with x0 = [r0; v0].
- `tspan::Tuple{TT,TT}`: The span of the independant variable forwhich to solve the ode.
- `mu::Real`: The mass parameter of the three-body system (`mu = m2 / (m1 + m2)`).

# keywords
- `solver::OrdinaryDiffEq.AbstractODESolver`: The solver to use for the integration.
- `reltol`: The relative tolerance for the solver.
- `abstol`: The absolute tolerance for the solver.

# returns
- `SMatrix{6,6,TX,36}`: The stm for the full trajectory
"""
@inline function propagate_return_final_state_and_stm(
    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
) where {TX,TT}
    # Form initial state with STM
    z_start = initial_state_with_stm(x0)

    # Propagate
    z_stop = solve(
        ODEProblem{false, SciMLBase.FullSpecialize}(
            natural_crtbp_eom_with_stm,
            z_start,
            tspan,
            (mu,),
        ),
        solver;
        abstol = abstol,
        reltol = reltol,
        save_start = false,
        save_everystep = false,
    )[end]

    return get_state_and_stm(z_stop)
end
