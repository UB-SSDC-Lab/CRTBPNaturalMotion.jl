"""
    AbstractManifold

Abstract type for manifolds.
"""
abstract type AbstractManifold end

"""
    ParameterListType

A Union type of the acceptable types for a list of parameter values.
"""
ParameterListType = Union{StepRangeLen, AbstractVector}

"""
    InvariantManifold{PO <: AbstractPeriodicOrbit}

A struct for invariant manifolds of periodic orbits.

# Fields
- `orbit::PO`: The periodic orbit.
- `Δr::Float64`: The manifold perturbation size (change in position) along stable/unstable
    eigenvector direction.
"""
struct InvariantManifold{PO <: AbstractPeriodicOrbit} <: AbstractManifold
    # The periodic orbit
    orbit::PO

    # The manifold perturbation size (change in position)
    Δr::Float64

    function InvariantManifold(orbit::PO, Δ::AbstractFloat) where {PO <: AbstractPeriodicOrbit}
        Δ < 0 && throw(ArgumentError("Δ must be non-negative"))
        return new{PO}(orbit, Float64(Δ))
    end
end

"""
    get_stable_manifold_initial_condition(
        man::InvariantManifold, left_pert::Bool, τ1::AbstractFloat, PT::ParameterizationType;
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    )

Returns the initial condition for a stable manifold trajectory as parameterized by `τ1` ∈ [0,1].

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1::AbstractFloat`: The parameter for the initial condition on the periodic orbit.
- `PT::ParameterizationType`: The parameterization type for the `τ1` parameter.

# Keyword Arguments
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.

# Returns
- `x0m::SVector{6,Float64}`: The initial condition for the stable manifold trajectory.
"""
function get_stable_manifold_initial_condition(
    man::InvariantManifold, left_pert::Bool, τ1::AbstractFloat, PT::ParameterizationType;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Get the initial state and monodromy matrix at τ1 point
    x0, M = get_state_and_monodromy_matrix(
        man.orbit, τ1, PT;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )

    # Perform eigen decomposition
    vals, vecs = eigen(M)

    # Determine the index of stable eigen value
    idx = argmin(abs.(real(vals)))

    # Get unit vector in direction of stable eigenvector
    stable_vec = real(vecs[:, idx])
    stable_vec /= norm(stable_vec)

    # If left_pert == true, we want to make sure we're perturbing the position in the
    # negative x-direction (i.e., to the left of the original position if looking from top down).
    if (left_pert && stable_vec[1] > 0.0) || (!left_pert && stable_vec[1] < 0.0)
        stable_vec = -stable_vec
    end

    # Compute the perturbation size that changes the position by Δr
    pert_sf = man.Δr / norm(stable_vec[SA[1,2,3]])

    # Return
    x0m = x0 + pert_sf*stable_vec
    return x0m
end

"""
    get_unstable_manifold_initial_condition(
        man::InvariantManifold, left_pert::Bool, τ1::AbstractFloat, PT::ParameterizationType;
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    )

Returns the initial condition for a unstable manifold trajectory as parameterized by `τ1` ∈ [0,1].

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1::AbstractFloat`: The parameter for the initial condition on the periodic orbit.
- `PT::ParameterizationType`: The parameterization type for the `τ1` parameter.

# Keyword Arguments
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.

# Returns
- `x0m::SVector{6,Float64}`: The initial condition for the stable manifold trajectory.
"""
function get_unstable_manifold_initial_condition(
    man::InvariantManifold, left_pert::Bool, τ1::AbstractFloat, PT::ParameterizationType;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Get the initial state and monodromy matrix at τ1 point
    x0, M = get_state_and_monodromy_matrix(
        man.orbit, τ1, PT;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )

    # Perform eigen decomposition
    vals, vecs = eigen(M)

    # Determine the index of unstable eigen value
    idx = argmax(abs.(real(vals)))

    # Get unit vector in direction of unstable eigenvector
    stable_vec = real(vecs[:, idx])
    stable_vec /= norm(stable_vec)

    # If left_pert == true, we want to make sure we're perturbing the position in the
    # negative x-direction (i.e., to the left of the original position if looking from top down).
    if (left_pert && stable_vec[1] > 0.0) || (!left_pert && stable_vec[1] < 0.0)
        stable_vec = -stable_vec
    end

    # Compute the perturbation size the position by Δr
    pert_sf = man.Δr / norm(stable_vec[SA[1,2,3]])

    # Return
    x0m = x0 + pert_sf*stable_vec
    return x0m
end

"""
    get_stable_manifold_trajectory(
        man::InvariantManifold, left_pert::Bool,
        τ1::AbstractFloat, PT1::ParameterizationType,
        τ2::Union{AbstractFloat,ParameterListType}, PT2::ParameterizationType;
        manifold_length::AbstractFloat = 1.0,
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    )

Returns a single stable manifold trajectory corresponding to the parameters τ1 and τ2.

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1::AbstractFloat`: The parameter for the initial condition on the periodic orbit.
- `PT1::ParameterizationType`: The parameterization type for the `τ1` parameter.
- `τ2::Union{AbstractFloat,ParameterListType}`: The parameter for the final condition on the
    manifold. If a single scalar is provided, the returned trajectory is a DifferentialEquations.jl
    solution struct. If a list is provided, the returned trajectory is an array of states
    corresponding to each element in the list.
- `PT2::ParameterizationType`: The parameterization type for the `τ2` parameter.

# Keyword Arguments
- `manifold_length::AbstractFloat = 1.0`: The length of the manifold, in CRTBP units. Sets
    the max length of the manifold, where `τ2 = 1.0` corresponds to the full `manifold_length`.
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.
"""
function get_stable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    τ2::AbstractFloat, PT2::ParameterizationType;
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Check that τ2 ∈ [0, 1]
    if τ2 < 0.0 || τ2 > 1.0
        error("τ2 must be in [0, 1]")
    end

    return propagate_return_all_states(
        get_stable_manifold_initial_condition(
            man, left_pert, τ1, PT1;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        (0.0, -τ2*manifold_length),
        man.orbit.mu,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end
function get_stable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    τ2::ParameterListType,  PT2::ParameterizationType;
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Check that τ2 ∈ [0, 1]
    if minimum(τ2) < 0.0 || maximum(τ2) > 1.0
        error("τ2 must be in [0, 1]")
    end

    return propagate_return_all_states(
        get_stable_manifold_initial_condition(
            man, left_pert, τ1, PT1;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        -τ2*manifold_length,
        man.orbit.mu,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

"""
    get_stable_manifold_trajectory(
        man::InvariantManifold, left_pert::Bool,
        τ1::AbstractFloat, PT1::ParameterizationType,
        start_cond::Union{AbstractFloat,Function},
        τ2::Union{AbstractFloat,ParameterListType}, PT2::ParameterizationType;
        manifold_length::AbstractFloat = 1.0,
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    )

Returns a single stable manifold trajectory corresponding to the parameters τ1 and τ2. Here,
the starting point of the returned trajectory is determined by the argument `start_cond`. If
`start_cond` is a scalar, this specifies the distance from the initial condition on the
periodic orbit in terms of `PT2`. If `start_cond` is a function, this specifies a termination
condition callback function for the solver that will be used to propagate to the initial state
of the returned trajectory from the initial condition on the periodic orbit.

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1::AbstractFloat`: The parameter for the initial condition on the periodic orbit.
- `PT1::ParameterizationType`: The parameterization type for the `τ1` parameter.
- `start_cond::Union{AbstractFloat,Function}`: The starting condition for the manifold trajectory.
    If a scalar, this specifies the distance from the initial condition on the periodic orbit in terms
    of `PT2`. If a function, this specifies a termination condition callback function for the solver
    that will be used to propagate to the initial state of the returned trajectory from the initial
    condition on the periodic orbit.
- `τ2::Union{AbstractFloat,ParameterListType}`: The parameter for the final condition on the
    manifold. If a single scalar is provided, the returned trajectory is a DifferentialEquations.jl
    solution struct. If a list is provided, the returned trajectory is an array of states
    corresponding to each element in the list.
- `PT2::ParameterizationType`: The parameterization type for the `τ2` parameter.

# Keyword Arguments
- `manifold_length::AbstractFloat = 1.0`: The length of the manifold, in CRTBP units. Sets
    the max length of the manifold, where `τ2 = 1.0` corresponds to the full `manifold_length`,
    starting from the state parameterized by `start_cond`.
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.
"""
function get_stable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    start_distance::AbstractFloat, τ2::AbstractFloat, PT2::ParameterizationType;
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    return propagate_return_all_states(
        propagate_return_final_state(
            get_stable_manifold_initial_condition(
                man, left_pert, τ1, PT1;
                solver = solver,
                reltol = reltol,
                abstol = abstol,
            ),
            (0.0, -start_distance),
            man.orbit.mu,
            PT2;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        (0.0, -τ2*manifold_length),
        man.orbit.mu,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = 1e-14,
    )
end
function get_stable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    start_distance::AbstractFloat, τ2::ParameterListType, PT2::ParameterizationType;
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    return propagate_return_all_states(
        propagate_return_final_state(
            get_stable_manifold_initial_condition(
                man, left_pert, τ1, PT1;
                solver = solver,
                reltol = reltol,
                abstol = abstol,
            ),
            (0.0, -start_distance),
            man.orbit.mu,
            PT2;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        -τ2*manifold_length,
        man.orbit.mu,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end
function get_stable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    start_cond::Function, τ2::AbstractFloat, PT2::ParameterizationType;
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    return propagate_return_all_states(
        propagate_return_final_state(
            get_stable_manifold_initial_condition(
                man, left_pert, τ1, PT1;
                solver = solver,
                reltol = reltol,
                abstol = abstol,
            ),
            (0.0, -Inf),
            man.orbit.mu,
            start_cond,
            PT2;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        (0.0, -τ2*manifold_length),
        man.orbit.mu,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end
function get_stable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    start_cond::Function, τ2::ParameterListType, PT2::ParameterizationType;
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    return propagate_return_all_states(
        propagate_return_final_state(
            get_stable_manifold_initial_condition(
                man, left_pert, τ1, PT1;
                solver = solver,
                reltol = reltol,
                abstol = abstol,
            ),
            (0.0, -Inf),
            man.orbit.mu,
            start_cond,
            PT2;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        -τ2*manifold_length,
        man.orbit.mu,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

"""
    get_stable_manifold_trajectory(
        man::InvariantManifold, left_pert::Bool,
        τ1::AbstractFloat,   PT1::ParameterizationType,
        term_cond::Function, PT2::ParameterizationType;
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    )

Returns a single stable manifold trajectory given scalar parameter τ1, where the final
state of the returned trajectory is determined through the satisfaction of `term_cond`, a
termination condition callback function.

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1::AbstractFloat`: The parameter for the initial condition on the periodic orbit.
- `PT1::ParameterizationType`: The parameterization type for the `τ1` parameter.
- `term_cond::Function`: The termination condition callback function for the solver.
- `PT2::ParameterizationType`: The parameterization type for the final condition on the manifold.

# Keyword Arguments
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.
"""
function get_stable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat,   PT1::ParameterizationType,
    term_cond::Function, PT2::ParameterizationType;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    return propagate_return_all_states(
        get_stable_manifold_initial_condition(
            man, left_pert, τ1, PT1;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        (0.0, -Inf),
        man.orbit.mu,
        term_cond,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

"""
    get_stable_manifold_trajectory(
        man::InvariantManifold, left_pert::Bool,
        τ1::AbstractFloat, PT1::ParameterizationType,
        start_cond::Function, term_cond::Function, PT2::ParameterizationType;
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    )

Returns a single stable manifold trajectory given scalar parameter τ1, where the starting
state of the returned trajectory is determined through the satisfaction of `start_cond` and
the final state is determined through the satisfaction of `term_cond`, where `start_cond` and
`term_cond` are termination condition callback functions.

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1::AbstractFloat`: The parameter for the initial condition on the periodic orbit.
- `PT1::ParameterizationType`: The parameterization type for the `τ1` parameter.
- `start_cond::Function`: The starting condition callback function for the solver.
- `term_cond::Function`: The termination condition callback function for the solver.
- `PT2::ParameterizationType`: The parameterization type for the final condition on the manifold.

# Keyword Arguments
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.
"""
function get_stable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    start_cond::Function, term_cond::Function, PT2::ParameterizationType;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    return propagate_return_all_states(
        propagate_return_final_state(
            get_stable_manifold_initial_condition(
                man, left_pert, τ1, PT1;
                solver = solver,
                reltol = reltol,
                abstol = abstol,
            ),
            (0.0, -Inf),
            man.orbit.mu,
            start_cond,
            PT2;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        (0.0, -Inf),
        man.orbit.mu,
        term_cond,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

"""
    get_unstable_manifold_trajectory(
        man::InvariantManifold, left_pert::Bool,
        τ1::AbstractFloat, PT1::ParameterizationType,
        τ2::Union{AbstractFloat,ParameterListType}, PT2::ParameterizationType;
        manifold_length::AbstractFloat = 1.0,
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    )

Returns a single unstable manifold trajectory corresponding to the parameters τ1 and τ2.

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1::AbstractFloat`: The parameter for the initial condition on the periodic orbit.
- `PT1::ParameterizationType`: The parameterization type for the `τ1` parameter.
- `τ2::Union{AbstractFloat,ParameterListType}`: The parameter for the final condition on the
    manifold. If a single scalar is provided, the returned trajectory is a DifferentialEquations.jl
    solution struct. If a list is provided, the returned trajectory is an array of states
    corresponding to each element in the list.
- `PT2::ParameterizationType`: The parameterization type for the `τ2` parameter.

# Keyword Arguments
- `manifold_length::AbstractFloat = 1.0`: The length of the manifold, in CRTBP units. Sets
    the max length of the manifold, where `τ2 = 1.0` corresponds to the full `manifold_length`.
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.
"""
function get_unstable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    τ2::AbstractFloat, PT2::ParameterizationType;
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Check that τ2 ∈ [0, 1]
    if τ2 < 0.0 || τ2 > 1.0
        error("τ2 must be in [0, 1]")
    end

    return propagate_return_all_states(
        get_unstable_manifold_initial_condition(
            man, left_pert, τ1, PT1;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        (0.0, τ2*manifold_length),
        man.orbit.mu,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end
function get_unstable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat,     PT1::ParameterizationType,
    τ2::ParameterListType, PT2::ParameterizationType;
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Check that τ2 ∈ [0, 1]
    if τ2 < 0.0 || τ2 > 1.0
        error("τ2 must be in [0, 1]")
    end

    return propagate_return_all_states(
        get_unstable_manifold_initial_condition(
            man, left_pert, τ1, PT1;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        τ2*manifold_length,
        man.orbit.mu,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

"""
    get_unstable_manifold_trajectory(
        man::InvariantManifold, left_pert::Bool,
        τ1::AbstractFloat, PT1::ParameterizationType,
        start_cond::Union{AbstractFloat,Function},
        τ2::Union{AbstractFloat,ParameterListType}, PT2::ParameterizationType;
        manifold_length::AbstractFloat = 1.0,
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    )

Returns a single unstable manifold trajectory corresponding to the parameters τ1 and τ2. Here,
the starting point of the returned trajectory is determined by the argument `start_cond`. If
`start_cond` is a scalar, this specifies the distance from the initial condition on the
periodic orbit in terms of `PT2`. If `start_cond` is a function, this specifies a termination
condition callback function for the solver that will be used to propagate to the initial state
of the returned trajectory from the initial condition on the periodic orbit.

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1::AbstractFloat`: The parameter for the initial condition on the periodic orbit.
- `PT1::ParameterizationType`: The parameterization type for the `τ1` parameter.
- `start_cond::Union{AbstractFloat,Function}`: The starting condition for the manifold trajectory.
    If a scalar, this specifies the distance from the initial condition on the periodic orbit in terms
    of `PT2`. If a function, this specifies a termination condition callback function for the solver
    that will be used to propagate to the initial state of the returned trajectory from the initial
    condition on the periodic orbit.
- `τ2::Union{AbstractFloat,ParameterListType}`: The parameter for the final condition on the
    manifold. If a single scalar is provided, the returned trajectory is a DifferentialEquations.jl
    solution struct. If a list is provided, the returned trajectory is an array of states
    corresponding to each element in the list.
- `PT2::ParameterizationType`: The parameterization type for the `τ2` parameter.

# Keyword Arguments
- `manifold_length::AbstractFloat = 1.0`: The length of the manifold, in CRTBP units. Sets
    the max length of the manifold, where `τ2 = 1.0` corresponds to the full `manifold_length`,
    starting from the state parameterized by `start_cond`.
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.
"""
function get_unstable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    start_distance::AbstractFloat, τ2::AbstractFloat, PT2::ParameterizationType;
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Propagate and return
    return propagate_return_all_states(
        propagate_return_final_state(
            get_unstable_manifold_initial_condition(
                man, left_pert, τ1, PT1;
                solver = solver,
                reltol = reltol,
                abstol = abstol,
            ),
            (0.0, start_distance),
            man.orbit.mu,
            PT2,
        ),
        (0.0, τ2*manifold_length),
        man.orbit.mu,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end
function get_unstable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    start_distance::AbstractFloat, τ2::ParameterListType, PT2::ParameterizationType;
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Propagate and return
    return propagate_return_all_states(
        propagate_return_final_state(
            get_unstable_manifold_initial_condition(
                man, left_pert, τ1, PT1;
                solver = solver,
                reltol = reltol,
                abstol = abstol,
            ),
            (0.0, start_distance),
            man.orbit.mu,
            PT2;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        τ2*manifold_length,
        man.orbit.mu,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end
function get_unstable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    start_cond::Function, τ2::AbstractFloat, PT2::ParameterizationType;
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Propagate and return
    return propagate_return_all_states(
        propagate_return_final_state(
            get_unstable_manifold_initial_condition(
                man, left_pert, τ1, PT1;
                solver = solver,
                reltol = reltol,
                abstol = abstol,
            ),
            (0.0, Inf),
            man.orbit.mu,
            start_cond,
            PT2;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        (0.0, τ2*manifold_length),
        man.orbit.mu,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end
function get_unstable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    start_cond::Function, τ2::ParameterListType, PT2::ParameterizationType;
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Propagate and return
    return propagate_return_all_states(
        propagate_return_final_state(
            get_unstable_manifold_initial_condition(
                man, left_pert, τ1, PT1;
                solver = solver,
                reltol = reltol,
                abstol = abstol,
            ),
            (0.0, Inf),
            man.orbit.mu,
            start_cond,
            PT2;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        τ2*manifold_length,
        man.orbit.mu,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

"""
    get_unstable_manifold_trajectory(
        man::InvariantManifold, left_pert::Bool,
        τ1::AbstractFloat,   PT1::ParameterizationType,
        term_cond::Function, PT2::ParameterizationType;
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    )

Returns a single unstable manifold trajectory given scalar parameter τ1, where the final
state of the returned trajectory is determined through the satisfaction of `term_cond`, a
termination condition callback function.

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1::AbstractFloat`: The parameter for the initial condition on the periodic orbit.
- `PT1::ParameterizationType`: The parameterization type for the `τ1` parameter.
- `term_cond::Function`: The termination condition callback function for the solver.
- `PT2::ParameterizationType`: The parameterization type for the final condition on the manifold.

# Keyword Arguments
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.
"""
function get_unstable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat,   PT1::ParameterizationType,
    term_cond::Function, PT2::ParameterizationType;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    return propagate_return_all_states(
        get_unstable_manifold_initial_condition(
            man, left_pert, τ1, PT1;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        (0.0, Inf),
        man.orbit.mu,
        term_cond,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

"""
    get_unstable_manifold_trajectory(
        man::InvariantManifold, left_pert::Bool,
        τ1::AbstractFloat, PT1::ParameterizationType,
        start_cond::Function, term_cond::Function, PT2::ParameterizationType;
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    )

Returns a single unstable manifold trajectory given scalar parameter τ1, where the starting
state of the returned trajectory is determined through the satisfaction of `start_cond` and
the final state is determined through the satisfaction of `term_cond`, where `start_cond` and
`term_cond` are termination condition callback functions.

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1::AbstractFloat`: The parameter for the initial condition on the periodic orbit.
- `PT1::ParameterizationType`: The parameterization type for the `τ1` parameter.
- `start_cond::Function`: The starting condition callback function for the solver.
- `term_cond::Function`: The termination condition callback function for the solver.
- `PT2::ParameterizationType`: The parameterization type for the final condition on the manifold.

# Keyword Arguments
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.
"""
function get_unstable_manifold_trajectory(
    man::InvariantManifold, left_pert::Bool,
    τ1::AbstractFloat, PT1::ParameterizationType,
    start_cond::Function, term_cond::Function, PT2::ParameterizationType;
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    return propagate_return_all_states(
        propagate_return_final_state(
            get_unstable_manifold_initial_condition(
                man, left_pert, τ1, PT1;
                solver = solver,
                reltol = reltol,
                abstol = abstol,
            ),
            (0.0, Inf),
            man.orbit.mu,
            start_cond,
            PT2;
            solver = solver,
            reltol = reltol,
            abstol = abstol,
        ),
        (0.0, Inf),
        man.orbit.mu,
        term_cond,
        PT2;
        solver = solver,
        reltol = reltol,
        abstol = abstol,
    )
end

"""
    generate_stable_manifold_cross_section_cheb_interpolant(
        man::InvariantManifold, left_pert::Bool,
        τ1_order::Int, PT1::ParameterizationType,
        τ2::AbstractFloat, PT2::ParameterizationType;
        start_cond::Union{Nothing,Function,AbstractFloat} = nothing,
        manifold_length::AbstractFloat = 1.0,
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) -> FastChebInterpolation{FastChebInterp.ChebPoly{1, SVector{6, Float64}, Int64}}

Generate a Chebyshev interpolant for a cross-section of the stable `InvariantManifold`.

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1_order::Int`: The order of the Chebyshev interpolation in the `τ1` direction.
- `PT1::ParameterizationType`: The parameterization type for the `τ1` parameter.
- `τ2::AbstractFloat`: The fixed choice for the `τ2` parameter, which specified the location
    of the fixed manifold cross-section.
- `PT2::ParameterizationType`: The parameterization type for the `τ2` parameter.

# Keyword Arguments
- `start_cond::Union{Nothing,Function,AbstractFloat} = nothing`: The starting condition for the
    manifold trajectory. If a scalar, this specifies the distance from the initial condition on the
    periodic orbit in terms of `PT2`. If a function, this specifies a termination condition callback
    function for the solver that will be used to propagate to the initial state of the returned
    trajectory from the initial condition on the periodic orbit.
- `manifold_length::AbstractFloat = 1.0`: The length of the manifold, in CRTBP units. Sets
    the max length of the manifold, where `τ2 = 1.0` corresponds to the full `manifold_length`.
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.

# Returns
- `FastChebInterpolation{FastChebInterp.ChebPoly{2, SVector{6, Float64}, Int64}}`: The Chebyshev
    interpolant for the stable manifold.
"""
function generate_stable_manifold_cross_section_cheb_interpolant(
    man::InvariantManifold, left_pert::Bool,
    τ1_order::Int, PT1::ParameterizationType,
    τ2::AbstractFloat, PT2::ParameterizationType;
    start_cond::Union{Nothing,Function,AbstractFloat} = nothing,
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Generate range of parameter values
    τ1s = chebpoints(τ1_order, 0.0, 1.0)

    # Loop over parameter values
    states = Vector{SVector{6,Float64}}(undef, τ1_order+1)
    for (i, τ1) in enumerate(τ1s)
        # Get manifold states
        states[i] = get_stable_manifold_trajectory(
            man, left_pert, τ1, PT1,
            start_cond, τ2, PT2;
            manifold_length = manifold_length,
            solver          = solver,
            reltol          = reltol,
            abstol          = abstol,
        )[end]
    end

    return FastChebInterpolation(chebinterp(states, 0.0, 1.0))
end

"""
    generate_stable_manifold_cross_section_cheb_approximation(
        man::InvariantManifold, left_pert::Bool,
        τ1_order::Int, PT1::ParameterizationType,
        τ2::AbstractFloat, PT2::ParameterizationType;
        start_cond::Union{Nothing,Function,AbstractFloat} = nothing,
        manifold_length::AbstractFloat = 1.0,
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) -> FastChebInterpolation{FastChebInterp.ChebPoly{1, SVector{6, Float64}, Int64}}

Generate a Chebyshev interpolant for a cross-section of the stable `InvariantManifold` through
a least-squares fit to

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1_npoints::Int`: The number of points to use in the `τ1` direction.
- `τ1_order::Int`: The order of the Chebyshev interpolation in the `τ1` direction.
- `PT1::ParameterizationType`: The parameterization type for the `τ1` parameter.
- `τ2::AbstractFloat`: The fixed choice for the `τ2` parameter, which specified the location
    of the fixed manifold cross-section.
- `PT2::ParameterizationType`: The parameterization type for the `τ2` parameter.

# Keyword Arguments
- `start_cond::Union{Nothing,Function,AbstractFloat} = nothing`: The starting condition for the
    manifold trajectory. If a scalar, this specifies the distance from the initial condition on the
    periodic orbit in terms of `PT2`. If a function, this specifies a termination condition callback
    function for the solver that will be used to propagate to the initial state of the returned
    trajectory from the initial condition on the periodic orbit.
- `manifold_length::AbstractFloat = 1.0`: The length of the manifold, in CRTBP units. Sets
    the max length of the manifold, where `τ2 = 1.0` corresponds to the full `manifold_length`.
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.

# Returns
- `FastChebInterpolation{FastChebInterp.ChebPoly{2, SVector{6, Float64}, Int64}}`: The Chebyshev
    interpolant for the stable manifold.
"""
function generate_stable_manifold_cross_section_cheb_approximation(
    man::InvariantManifold, left_pert::Bool,
    τ1_npoints::Int, τ1_order::Int, PT1::ParameterizationType,
    τ2::AbstractFloat, PT2::ParameterizationType;
    start_cond::Union{Nothing,Function,AbstractFloat} = nothing,
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Generate range of parameter values
    τ1s = chebpoints(τ1_npoints - 1, 0.0, 1.0)

    # Loop over parameter values
    states = Vector{SVector{6,Float64}}(undef, τ1_npoints)
    for (i, τ1) in enumerate(τ1s)
        # Get manifold states
        states[i] = get_stable_manifold_trajectory(
            man, left_pert, τ1, PT1,
            start_cond, τ2, PT2;
            manifold_length = manifold_length,
            solver          = solver,
            reltol          = reltol,
            abstol          = abstol,
        )[end]
    end

    return FastChebInterpolation(chebregression(τ1s, states, 0.0, 1.0, τ1_order))
end

"""
    generate_stable_manifold_cheb_interpolant(
        man::InvariantManifold, left_pert::Bool,
        τ1_order::Int, PT1::ParameterizationType,
        τ2_order::Int, PT2::ParameterizationType;
        start_cond::Union{Nothing,Function,AbstractFloat} = nothing,
        manifold_length::AbstractFloat = 1.0,
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) -> FastChebInterpolation{FastChebInterp.ChebPoly{2, SVector{6, Float64}, Int64}}

Generate a Chebyshev interpolant for the stable `InvariantManifold`.

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1_order::Int`: The order of the Chebyshev interpolation in the `τ1` direction.
- `PT1::ParameterizationType`: The parameterization type for the `τ1` parameter.
- `τ2_order::Int`: The order of the Chebyshev interpolation in the `τ2` direction.
- `PT2::ParameterizationType`: The parameterization type for the `τ2` parameter.

# Keyword Arguments
- `start_cond::Union{Nothing,Function,AbstractFloat} = nothing`: The starting condition for the
    manifold trajectory. If a scalar, this specifies the distance from the initial condition on the
    periodic orbit in terms of `PT2`. If a function, this specifies a termination condition callback
    function for the solver that will be used to propagate to the initial state of the returned
    trajectory from the initial condition on the periodic orbit.
- `manifold_length::AbstractFloat = 1.0`: The length of the manifold, in CRTBP units. Sets
    the max length of the manifold, where `τ2 = 1.0` corresponds to the full `manifold_length`.
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.

# Returns
- `FastChebInterpolation{FastChebInterp.ChebPoly{2, SVector{6, Float64}, Int64}}`: The Chebyshev
    interpolant for the stable manifold.
"""
function generate_stable_manifold_cheb_interpolant(
    man::InvariantManifold, left_pert::Bool,
    τ1_order::Int, PT1::ParameterizationType,
    τ2_order::Int, PT2::ParameterizationType;
    start_cond::Union{Nothing,Function,AbstractFloat} = nothing,
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Generate range of parameter values
    lb = SA[0,0]
    ub = SA[1,1]
    τs = chebpoints(SA[τ1_order, τ2_order], lb, ub)

    # These are avoidable allocations, but we only call this once
    # before using the interp a bunch
    τ1s = [τs[i,1][1] for i in axes(τs, 1)]
    τ2s = [τs[1,end - i + 1][2] for i in axes(τs, 2)]

    # Loop over parameter values
    states = Matrix{SVector{6,Float64}}(undef, τ1_order+1, τ2_order+1)
    for (i, τ1) in enumerate(τ1s)
        # Get manifold states
        τ1_traj = get_stable_manifold_trajectory(
            man, left_pert, τ1, PT1,
            start_cond, τ2s, PT2;
            manifold_length = manifold_length,
            solver          = solver,
            reltol          = reltol,
            abstol          = abstol,
        )
        for j in eachindex(τ2s)
            states[i,j] = τ1_traj[end - j + 1]
        end
    end

    return FastChebInterpolation(chebinterp(states, lb, ub))
end

"""
    generate_stable_manifold_cheb_approximation(
        man::InvariantManifold, left_pert::Bool,
        τ1_npoints::Int, τ1_order::Int, PT1::ParameterizationType,
        τ2_npoints::Int, τ2_order::Int, PT2::ParameterizationType;
        start_cond::Union{Nothing,Function,AbstractFloat} = nothing,
        manifold_length::AbstractFloat = 1.0,
        solver = Vern9(),
        reltol = 1e-14,
        abstol = 1e-14,
    ) -> FastChebInterpolation{FastChebInterp.ChebPoly{2, SVector{6, Float64}, Int64}}

Generate a Chebyshev approximation for the stable `InvariantManifold` in terms of `PT1`
and `PT2` for the parameters `τ1` and `τ2` employing `τ1_npoints` and `τ2_npoints` to a
Chebyshev polynomial of order `τ1_order` and `τ2_order`, respectively.

# Arguments
- `man::InvariantManifold`: The invariant manifold.
- `left_pert::Bool`: If true, the perturbation is to the left of the periodic orbit.
- `τ1_npoints::Int`: The number of points to use in the `τ1` direction.
- `τ1_order::Int`: The order of the Chebyshev interpolation in the `τ1` direction.
- `PT1::ParameterizationType`: The parameterization type for the `τ1` parameter.
- `τ2_npoints::Int`: The number of points to use in the `τ2` direction.
- `τ2_order::Int`: The order of the Chebyshev interpolation in the `τ2` direction.
- `PT2::ParameterizationType`: The parameterization type for the `τ2` parameter.

# Keyword Arguments
- `start_cond::Union{Nothing,Function,AbstractFloat} = nothing`: The starting condition for the
    manifold trajectory. If a scalar, this specifies the distance from the initial condition on the
    periodic orbit in terms of `PT2`. If a function, this specifies a termination condition callback
    function for the solver that will be used to propagate to the initial state of the returned
    trajectory from the initial condition on the periodic orbit.
- `manifold_length::AbstractFloat = 1.0`: The length of the manifold, in CRTBP units. Sets
    the max length of the manifold, where `τ2 = 1.0` corresponds to the full `manifold_length`.
- `solver = Vern9()`: The DifferentialEquations.jl solver to use.
- `reltol = 1e-14`: The relative tolerance for the solver.
- `abstol = 1e-14`: The absolute tolerance for the solver.

# Returns
- `FastChebInterpolation{FastChebInterp.ChebPoly{2, SVector{6, Float64}, Int64}}`: The Chebyshev
    interpolant for the stable manifold.
"""
function generate_stable_manifold_cheb_approximation(
    man::InvariantManifold, left_pert::Bool,
    τ1_npoints::Int, τ1_order::Int, PT1::ParameterizationType,
    τ2_npoints::Int, τ2_order::Int, PT2::ParameterizationType;
    start_cond::Union{Nothing,Function,AbstractFloat} = nothing,
    manifold_length::AbstractFloat = 1.0,
    solver = Vern9(),
    reltol = 1e-14,
    abstol = 1e-14,
)
    # Generate range of parameter values
    lb = SA[0,0]
    ub = SA[1,1]
    τs = chebpoints(SA[τ1_npoints - 1, τ2_npoints - 1], lb, ub)

    # These are avoidable allocations, but we only call this once
    # before using the interp a bunch
    τ1s = [τs[i,1][1] for i in axes(τs, 1)]
    τ2s = [τs[1,end - i + 1][2] for i in axes(τs, 2)]

    # Loop over parameter values
    idx     = 0
    xs      = Vector{SVector{2,Float64}}(undef, τ1_npoints*τ2_npoints)
    states  = Vector{SVector{6,Float64}}(undef, τ1_npoints*τ2_npoints)
    for (i, τ1) in enumerate(τ1s)
        # Get manifold states
        τ1_traj = get_stable_manifold_trajectory(
            man, left_pert, τ1, PT1,
            start_cond, τ2s, PT2;
            manifold_length = manifold_length,
            solver          = solver,
            reltol          = reltol,
            abstol          = abstol,
        )
        for j in eachindex(τ2s)
            idx += 1
            xs[idx] = SA[τ1, τ2s[j]]
            states[idx] = τ1_traj[end - j + 1]
        end
    end

    return FastChebInterpolation(chebregression(xs, states, lb, ub, (τ1_order, τ2_order)))
end
