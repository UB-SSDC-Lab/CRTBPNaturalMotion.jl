var documenterSearchIndex = {"docs":
[{"location":"lib/internal/temp/#Internals","page":"Internals","title":"Internals","text":"","category":"section"},{"location":"lib/internal/temp/","page":"Internals","title":"Internals","text":"Listing all non-exported types and functions here for now, but split off categories onto separate pages in the future!","category":"page"},{"location":"lib/internal/temp/","page":"Internals","title":"Internals","text":"Modules = [CRTBPNaturalMotion]\nPublic = false","category":"page"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.ParameterizationType","page":"Internals","title":"CRTBPNaturalMotion.ParameterizationType","text":"ParameterizationType\n\nA Union type for the different types of parameterizations that can be used when parameterizing a periodic orbit or its associated invairiant stable or intable manifold.\n\n\n\n\n\n","category":"type"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.AbstractIndependantVariable","page":"Internals","title":"CRTBPNaturalMotion.AbstractIndependantVariable","text":"AbstractIndependantVariable\n\nAn abstract type to encompas all supported independant variables.\n\n\n\n\n\n","category":"type"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.get_full_initial_state-Tuple{TypeAPeriodicOrbit}","page":"Internals","title":"CRTBPNaturalMotion.get_full_initial_state","text":"get_full_initial_state(orbit::TypeAPeriodicOrbit)\n\nGet the full initial reference state of the periodic orbit.\n\nArguments\n\norbit::TypeAPeriodicOrbit: The periodic orbit.\n\nReturns\n\nx0::SVector{6,Float64}: The full initial state of the periodic orbit.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.get_jacobi_integral-Tuple{CRTBPNaturalMotion.AbstractPeriodicOrbit}","page":"Internals","title":"CRTBPNaturalMotion.get_jacobi_integral","text":"get_jacobi_integral(orbit::AbstractPeriodicOrbit)\n\nGet the Jacobi integral of the periodic orbit.\n\nArguments\n\norbit::AbstractPeriodicOrbit: The periodic orbit.\n\nReturns\n\nC::Float64: The Jacobi integral of the periodic orbit.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.get_orbit_length-Tuple{CRTBPNaturalMotion.AbstractPeriodicOrbit, Type{Time}}","page":"Internals","title":"CRTBPNaturalMotion.get_orbit_length","text":"get_orbit_length(orbit::TypeAPeriodicOrbit, type::Type{AbstractIndependantVariable})\n\nGet the length of the periodic orbit in the desired variable.\n\nArguments\n\norbit::TypeAPeriodicOrbit: The periodic orbit.\ntype::Type{AbstractIndependantVariable}: The type of length to return.\n\nReturns\n\nlength::Float64: The length of the periodic orbit in the desired variable.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.get_state_and_monodromy_matrix-Tuple{CRTBPNaturalMotion.AbstractPeriodicOrbit, AbstractFloat, Union{Type{ArcLength}, Type{Time}}}","page":"Internals","title":"CRTBPNaturalMotion.get_state_and_monodromy_matrix","text":"get_state_and_monodromy_matrix(orbit::AbstractPeriodicOrbit, τ1::AbstractFloat, PT::ParameterizationType)\n\nGet the state and monodromy matrix of the periodic orbit at a give state parameterized by τ1 ∈ [0, 1] in terms of PT.\n\nArguments\n\norbit::AbstractPeriodicOrbit: The periodic orbit.\nτ1::AbstractFloat: The parameterization of the final state on the orbit.\nPT::ParameterizationType: The type of parameterization to use.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.get_state_and_stm-Union{Tuple{SMatrix{6, 7, T, 42}}, Tuple{T}} where T","page":"Internals","title":"CRTBPNaturalMotion.get_state_and_stm","text":"get_state_and_stm(z::SMatrix{6,7,T,42}) where T\n\nGet the state and state transition matrix (STM) from the full integration matrix.\n\nArguments\n\nz::SMatrix{6,7,T,42}: Matrix with first column as the state and the remaining 6 x 6 as   the STM.\n\nReturns\n\nTuple{SVector{6,T},SMatrix{6,6,T,36}: The state vector and the state transition matrix   in a return tuple.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.get_stm-Union{Tuple{SMatrix{6, 7, T, 42}}, Tuple{T}} where T","page":"Internals","title":"CRTBPNaturalMotion.get_stm","text":"get_stm(z::SMatrix{6,7,T,42}) where T\n\nGet the state transition matrix (STM) from the full integration matrix.\n\nArguments\n\nz::SMatrix{6,7,T,42}: Matrix with first column as the state and the remaining 6 x 6 as   the STM.\n\nReturns\n\nSMatrix{6,6,T,36}: The state transition matrix.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.initial_state_with_stm-Union{Tuple{SVector{6, T}}, Tuple{T}} where T","page":"Internals","title":"CRTBPNaturalMotion.initial_state_with_stm","text":"initial_state_with_stm(x0::SVector{6,T}) where T\n\nCreate the initial state with the state transition matrix (STM).\n\nArguments\n\nx0::SVector{6,T}: Initial state vector.\n\nReturns\n\nSMatrix{6,7,T,42}: Matrix with first column as the state and the remaining 6 x 6 as   the STM initial condition (i.e., [x0, I]).\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.jacobi_integral_gradient-Tuple{Any, Any}","page":"Internals","title":"CRTBPNaturalMotion.jacobi_integral_gradient","text":"jacobi_integral_gradient(x, μ)\n\nComputes the gradient of the Jacobi integral for the CRTBP given the state x and mass parameter μ.\n\nArguments\n\nx::AbstractVector: State vector.\nμ::Real: Mass parameter.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.natural_crtbp_and_time_eom_with_independant_arclen-Tuple{Any, Any, Any}","page":"Internals","title":"CRTBPNaturalMotion.natural_crtbp_and_time_eom_with_independant_arclen","text":"natural_crtbp_and_time_eom_with_independant_arclen(x, p, t)\n\nCompute the equations of motion for the CRTBP model with an extra state  for the time with only natural motion. Here the independant variable has been  transformed from time to arc-length.\n\nArguments\n\nx::SVector{7}: State vector, i.e., [position^T; velocity^T].\np::Tuple{Float64,}: Tuple of parameters, (mu,)\nt::Float64: Time.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.natural_crtbp_eom-Tuple{Any, Any, Any}","page":"Internals","title":"CRTBPNaturalMotion.natural_crtbp_eom","text":"natural_crtbp_eom(x, p, t)\n\nCompute the equations of motion for the CRTBP model with only natural motion.\n\nArguments\n\nx::SVector{6}: State vector, i.e., [position^T; velocity^T].\np::Tuple{Float64,Float64,Float64}: Tuple of parameters, (mu,)\nt::Float64: Time.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.natural_crtbp_eom_with_arclen-Tuple{Any, Any, Any}","page":"Internals","title":"CRTBPNaturalMotion.natural_crtbp_eom_with_arclen","text":"natural_crtbp_eom_with_arclen(x, p, t)\n\nCompute the equations of motion for the CRTBP model with only natural motion,  employing the arclength along the trajectory as an additional state.\n\nArguments\n\nx::SVector{7}: State vector, i.e., [position^T; velocity^T; arc-len].\np::Tuple{Float64,}: Tuple of parameters, (mu,)\nt::Float64: Time.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.natural_crtbp_eom_with_independant_arclen-Tuple{Any, Any, Any}","page":"Internals","title":"CRTBPNaturalMotion.natural_crtbp_eom_with_independant_arclen","text":"natural_crtbp_eom_with_independant_arclen(x, p, t)\n\nCompute the equations of motion for the CRTBP model with only natural motion where the independant variable has been transformed from time to arc-length.\n\nArguments\n\nx::SVector{7}: State vector, i.e., [position^T; velocity^T].\np::Tuple{Float64,}: Tuple of parameters, (mu,)\nt::Float64: Time.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.natural_crtbp_eom_with_stm-Tuple{Any, Any, Any}","page":"Internals","title":"CRTBPNaturalMotion.natural_crtbp_eom_with_stm","text":"natural_crtbp_eom_with_stm(x, p, t)\n\nCompute the equations of motion for the CRTBP model with the state transition matrix.\n\nArguments\n\nx::SMatrix{6,7,Float64}: Matrix with first column ast the state and the remaining 6 x 6 as the STM.\np::Tuple{Float64,}: Tuple of parameters, (mu,)\nt::Float64: Time.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.propagate_from_initial_conditions-Tuple{CRTBPNaturalMotion.AbstractPeriodicOrbit, AbstractFloat, Union{Type{ArcLength}, Type{Time}}}","page":"Internals","title":"CRTBPNaturalMotion.propagate_from_initial_conditions","text":"propagate_from_initial_conditions(\n    orbit::AbstractPeriodicOrbit, τ1::AbstractFloat, PT::ParameterizationType,\n    solver = Vern9(),\n    reltol = 1e-14,\n    abstol = 1e-14,\n)\n\nPropagate the periodic orbit from the initial conditions to another state on the orbit, as parameterized by τ1 ∈ [0, 1] in terms of PT::ParameterizationType.\n\nArguments\n\norbit::AbstractPeriodicOrbit: The periodic orbit.\nτ1::AbstractFloat: The parameterization of the final state on the orbit.\nPT::ParameterizationType: The type of parameterization to use.\n\nKeyword Arguments\n\nsolver::OrdinaryDiffEq.AbstractODESolver: The ODE solver to use for the integration.\nreltol::Real: The relative tolerance for the ODE solver.\nabstol::Real: The absolute tolerance for the ODE solver.\n\nReturns\n\nx1::SVector{6,Float64}: The final state on the periodic orbit.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.propagate_return_final_state_and_stm-Union{Tuple{TT}, Tuple{TX}, Tuple{SVector{6, TX}, Tuple{TT, TT}, Any}} where {TX, TT}","page":"Internals","title":"CRTBPNaturalMotion.propagate_return_final_state_and_stm","text":"propagate_return_final_state_and_stm(\n    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu;\n    solver = Vern9(),\n    reltol = 1e-14,\n    abstol = 1e-14,\n) where {TX,TT}\n\nPropagate the state x0 over the time span tspan and return the final state and the state transition matrix (STM).\n\narguments\n\nx0::SVector{6,TX}: The initial state with x0 = [r0; v0].\ntspan::Tuple{TT,TT}: The span of the independant variable forwhich to solve the ode.\nmu::Real: The mass parameter of the three-body system (mu = m2 / (m1 + m2)).\n\nkeywords\n\nsolver::OrdinaryDiffEq.AbstractODESolver: The solver to use for the integration.\nreltol: The relative tolerance for the solver.\nabstol: The absolute tolerance for the solver.\n\nreturns\n\nSMatrix{6,6,TX,36}: The stm for the full trajectory\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.typeA_shooting_function-Union{Tuple{T}, Tuple{SVector{4, T}, Any}} where T","page":"Internals","title":"CRTBPNaturalMotion.typeA_shooting_function","text":"typeA_shooting_function(u, p; constraint = (:x_start_coordinate, 0.0))\n\nShooting function for computing Type A periodic orbits in the CRTBP using the methodology described by Restrepo and Russell in A database of planar axisymmetric periodic orbits for the solar system doi: 10.1007/s10569-018-9844-6.\n\nArguments\n\nu::SVector{4,T}: Vector of decision variables, [rx, rz, vy, P / 2], where P is half of   the orbital period.\np::Tuple{Float64,}: Tuple of parameters, (mu,), where mu is the mass parameter.\n\nKeyword Arguments\n\nconstraint::Tuple{Symbol,Real}: Tuple of the constraint type and value corresponding   to the final constraint considered to provide a fully defined system of nonlinear   eqautions. The constraint type can be either :x_start_coordinate or :jacobi_integral.   The constraint value is the value that the constraint should be equal to for the desired   periodic orbit. Note that new constraints can be added with relative ease through   modification of this function and its corresponding jacobian function   typeA_shooting_jacobian.\nsolver::OrdinaryDiffEq.AbstractODESolver: The solver to use for the integration.\nreltol: The relative tolerance for the solver.\nabstol: The absolute tolerance for the solver.\n\nReturns\n\nSVector{4,T}: Vector of residuals for the nonlinear system of equations.\n\n\n\n\n\n","category":"method"},{"location":"lib/internal/temp/#CRTBPNaturalMotion.typeA_shooting_jacobian-Union{Tuple{T}, Tuple{SVector{4, T}, Any}} where T","page":"Internals","title":"CRTBPNaturalMotion.typeA_shooting_jacobian","text":"typeA_shooting_jacobian(u, p; constraint = (:x_start_coordinate, 0.0))\n\nThe Jacobian of the shooting function for computing Type A periodic orbits in the CRTBP using the methodology described by Restrepo and Russell in A database of planar axisymmetric periodic orbits for the solar system doi: 10.1007/s10569-018-9844-6.\n\nArguments\n\nu::SVector{4,T}: Vector of decision variables, [rx, rz, vy, P / 2], where P is half of   the orbital period.\np::Tuple{Float64,}: Tuple of parameters, (mu,), where mu is the mass parameter.\n\nKeyword Arguments\n\nconstraint::Tuple{Symbol,Real}: Tuple of the constraint type and value corresponding   to the final constraint considered to provide a fully defined system of nonlinear   eqautions. The constraint type can be either :x_start_coordinate or :jacobi_integral.\nsolver::OrdinaryDiffEq.AbstractODESolver: The solver to use for the integration.\nreltol: The relative tolerance for the solver.\nabstol: The absolute tolerance for the solver.\n\nReturns\n\nSMatrix{4,4,T,16}: Jacobian of the nonlinear system of equations.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#Public-API-Documentation","page":"Public API","title":"Public API Documentation","text":"","category":"section"},{"location":"lib/public/","page":"Public API","title":"Public API","text":"Documentation for CRTBPNaturalMotion's public interface.","category":"page"},{"location":"lib/public/#Contents","page":"Public API","title":"Contents","text":"","category":"section"},{"location":"lib/public/","page":"Public API","title":"Public API","text":"Pages = [\"public.md\"]\nDepth = 2:2","category":"page"},{"location":"lib/public/#Index","page":"Public API","title":"Index","text":"","category":"section"},{"location":"lib/public/","page":"Public API","title":"Public API","text":"Pages = [\"public.md\"]","category":"page"},{"location":"lib/public/#Periodic-Orbit-Computation","page":"Public API","title":"Periodic Orbit Computation","text":"","category":"section"},{"location":"lib/public/","page":"Public API","title":"Public API","text":"Modules = [CRTBPNaturalMotion]\nPages = [\"PeriodicOrbitStructs.jl\", \"periodic_orbit_correction.jl\"]\nPrivate = false","category":"page"},{"location":"lib/public/#CRTBPNaturalMotion.TypeAPeriodicOrbit","page":"Public API","title":"CRTBPNaturalMotion.TypeAPeriodicOrbit","text":"TypeAPeriodicOrbit\n\nA representation of a CRTBP periodic orbit of Type A.\n\nFields\n\nu0::SVector{3, Float64}: The initial conditions for the periodic orbit.\nP::Float64: The period of the periodic orbit.\nS::Float64: The arc length along the periodic orbit.\nN_cross::Int: The number of crossings through x-z plane in half period.\nmu::Float64: The gravitational parameter of the CRTBP.\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#CRTBPNaturalMotion.TypeAPeriodicOrbit-NTuple{4, Any}","page":"Public API","title":"CRTBPNaturalMotion.TypeAPeriodicOrbit","text":"TypeAPeriodicOrbit(\n    rx, rz, vy, mu;\n    N_cross = 1,\n    constraint = (:x_start_coordinate, 0.8),\n    ode_solver = Vern9(),\n    ode_reltol = 1e-12,\n    ode_abstol = 1e-12,\n    nl_solver  = SimpleTrustRegion(autodiff = nothing, nlsolve_update_rule = Val(true)),\n)\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#CRTBPNaturalMotion.TypeAPeriodicOrbit-Tuple{TypeAPeriodicOrbit}","page":"Public API","title":"CRTBPNaturalMotion.TypeAPeriodicOrbit","text":"TypeAPeriodicOrbit(\n    orbit::TypeAPeriodicOrbit;\n    constraint = (:x_start_coordinate, 0.8),\n    ode_solver = Vern9(),\n    ode_reltol = 1e-12,\n    ode_abstol = 1e-12,\n    nl_solver  = SimpleTrustRegion(autodiff = nothing, nlsolve_update_rule = Val(true)),\n)\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#CRTBPNaturalMotion.get_full_orbit-Tuple{CRTBPNaturalMotion.AbstractPeriodicOrbit}","page":"Public API","title":"CRTBPNaturalMotion.get_full_orbit","text":"get_full_orbit(orbit::AbstractPeriodicOrbit)\n\nGet the full orbit of the periodic orbit.\n\nArguments\n\norbit::AbstractPeriodicOrbit: The periodic orbit.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#CRTBPNaturalMotion.correct_typeA_initial_conditions-NTuple{4, Any}","page":"Public API","title":"CRTBPNaturalMotion.correct_typeA_initial_conditions","text":"correct_typeA_initial_conditions(\n    rx_guess, rz_guess, vy_guess, mu;\n    N_cross     = 1,\n    constraint  = (:x_start_coordinate, 0.0),\n    ode_solver  = Vern9(),\n    ode_reltol  = 1e-14,\n    ode_abstol  = 1e-14,\n    nl_solver   = SimpleTrustRegion(autodiff = nothing, nlsolve_update_rule = Val(true)),\n)\n\nCorrect the initial conditions of a type-A periodic orbit using a forward shooting method.\n\nArguments\n\nrx_guess::Real: Initial guess for the x-coordinate of the initial state.\nrz_guess::Real: Initial guess for the z-coordinate of the initial state.\nvy_guess::Real: Initial guess for the y-velocity of the initial state.\nmu::Real: Gravitational parameter of the CRTBP.\n\nKeyword Arguments\n\nN_cross::Int = 1: Number of crossings of the xz-plane to consider in the half period.\nconstraint::Tuple{Symbol,Real}: Tuple of the constraint type and value corresponding   to the final constraint considered to provide a fully defined system of nonlinear   eqautions. The constraint type can be either :x_start_coordinate or :jacobi_integral.   The constraint value is the value that the constraint should be equal to for the desired   periodic orbit.\node_solver::OrdinaryDiffEq.AbstractODESolver: The ODE solver to use for the integration.\node_reltol::Real: The relative tolerance for the ODE solver.\node_abstol::Real: The absolute tolerance for the ODE solver.\nnl_solver::NonlinearSolve.AbstractNLsolveSolver: The nonlinear solver to use for the   root-finding problem. Note, there seems to be a strange bug with the NonlinearSolve.jl   TrustRegion and NewtonRaphson solvers when using out-of-place Jacobian defined   to return an SMatrix. Therefore, it is recommened to use either the default   SimpleTrustRegion solver or the SimpleNewtonRaphson solver (which are specifically   tailored towards StaticArrays.jl arrays).\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#Type-Flags","page":"Public API","title":"Type Flags","text":"","category":"section"},{"location":"lib/public/","page":"Public API","title":"Public API","text":"Modules = [CRTBPNaturalMotion]\nPages = [\"type_flags.jl\"]\nPrivate = false","category":"page"},{"location":"lib/public/#CRTBPNaturalMotion.ArcLength","page":"Public API","title":"CRTBPNaturalMotion.ArcLength","text":"ArcLength\n\nA concrete type to represent arc-length as the independant variable.\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#CRTBPNaturalMotion.Time","page":"Public API","title":"CRTBPNaturalMotion.Time","text":"Time\n\nA concrete type to represent time as the independant variable.\n\n\n\n\n\n","category":"type"},{"location":"lib/public/#Propagation-Functions","page":"Public API","title":"Propagation Functions","text":"","category":"section"},{"location":"lib/public/","page":"Public API","title":"Public API","text":"Modules = [CRTBPNaturalMotion]\nPages = [\"integration.jl\"]\nPrivate = false","category":"page"},{"location":"lib/public/#CRTBPNaturalMotion.propagate_return_all_states-Union{Tuple{TT}, Tuple{TX}, Tuple{SVector{6, TX}, Tuple{TT, TT}, Any, Function, Function, Type{Time}}} where {TX, TT}","page":"Public API","title":"CRTBPNaturalMotion.propagate_return_all_states","text":"propagate_return_all_states(\n    x0::SVector{6,TX},\n    tspan::Tuple{TT,TT},\n    mu,\n    cond::Function,\n    affect::Function,\n    iv::Type{AbstractIndependantVariable};\n    solver = Vern9(),\n    reltol = 1e-14,\n    abstol = 1e-14,\n) where {TX,TT}\n\nPropagate the state x0 over the time span tspan and return the full ODE solution. This is done employing a DifferentialEquations.jl ContinuousCallback defined with cond and affect. The final argument can be specified to set the desired independant variable (i.e., Time or ArcLength).\n\nArguments\n\nx0::SVector{6,TX}: The initial state with x0 = [r0; v0].\ntspan::Tuple{TT,TT}: The span of the independant variable forwhich to solve the ODE.\nmu::Real: The mass parameter of the three-body system (mu = m2 / (m1 + m2)).\ncond::Function: The DifferentialEquations.jl continuous callback condition function   that will apply affect to the integrator struct (see   DifferentialEquations.jl callback documentation).\naffect::Function: The DifferentialEquations.jl continuous callback affect function.\niv::Type{AbstractIndependantVariable}: The desired independant variable (Time if   unspecified).\n\nKeywords\n\nsolver::OrdinaryDiffEq.AbstractODESolver: The solver to use for the integration.\nreltol: The relative tolerance for the solver.\nabstol: The absolute tolerance for the solver.\n\nReturns\n\nODESolution: The full solution to the ODE problem.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#CRTBPNaturalMotion.propagate_return_all_states-Union{Tuple{TT}, Tuple{TX}, Tuple{SVector{6, TX}, Tuple{TT, TT}, Any, Function, Type{Time}}} where {TX, TT}","page":"Public API","title":"CRTBPNaturalMotion.propagate_return_all_states","text":"propagate_return_all_states(\n    x0::SVector{6,TX},\n    tspan::Tuple{TT,TT},\n    mu,\n    term_cond::Function[, iv::Type{AbstractIndependantVariable}];\n    solver = Vern9(),\n    reltol = 1e-14,\n    abstol = 1e-14,\n) where {TX,TT}\n\nPropagate the state x0 over the time span tspan and return the full ODE solution. At the first point at which term_cond is zero (if any), the integration will terminate. The final argument can be specified to set the desired independant variable (i.e., Time or ArcLength).\n\nArguments\n\nx0::SVector{6,TX}: The initial state with x0 = [r0; v0].\ntspan::Tuple{TT,TT}: The span of the independant variable forwhich to solve the ODE.\nmu::Real: The mass parameter of the three-body system (mu = m2 / (m1 + m2)).\nterm_cond::Function: The DifferentialEquations.jl continuous callback condition function   that will terminate the integration when it is zero. Function should be of the form   term_cond(x,t,integ), where x is the state, t is the independant variable, and   integ is the integrator struct (see DifferentialEquations.jl callback documentation).\niv::Type{AbstractIndependantVariable}: The desired independant variable (Time if   unspecified).\n\nKeywords\n\nsolver::OrdinaryDiffEq.AbstractODESolver: The solver to use for the integration.\nreltol: The relative tolerance for the solver.\nabstol: The absolute tolerance for the solver.\n\nReturns\n\nODESolution: The full solution to the ODE problem.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#CRTBPNaturalMotion.propagate_return_all_states-Union{Tuple{TT}, Tuple{TX}, Tuple{SVector{6, TX}, Tuple{TT, TT}, Any, Type{Time}}} where {TX, TT}","page":"Public API","title":"CRTBPNaturalMotion.propagate_return_all_states","text":"propagate_return_all_states(\n    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu [, iv::Type{AbstractIndependantVariable}];\n    solver = Vern9(),\n    reltol = 1e-14,\n    abstol = 1e-14,\n) where {TX,TT}\n\nPropagate the state x0 over the time span tspan and return the full ODE solution. The final argument can be specified to set the desired independant variable (i.e., Time or ArcLength).\n\nArguments\n\nx0::SVector{6,T}: The initial state with x0 = [r0; v0].\ntspan::Tuple{T,T}: The span of the independant variable forwhich to solve the ODE.\nmu::Real: The mass parameter of the three-body system (mu = m2 / (m1 + m2)).\niv::Type{AbstractIndependantVariable}: The desired independant variable (Time if   unspecified).\n\nKeywords\n\nsolver::OrdinaryDiffEq.AbstractODESolver: The solver to use for the integration.\nreltol: The relative tolerance for the solver.\nabstol: The absolute tolerance for the solver.\n\nReturns\n\nODESolution: The full solution to the ODE problem.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#CRTBPNaturalMotion.propagate_return_all_states-Union{Tuple{TX}, Tuple{SVector{6, TX}, AbstractVector, Any, Type{Time}}} where TX","page":"Public API","title":"CRTBPNaturalMotion.propagate_return_all_states","text":"propagate_return_all_states(\n    x0::SVector{6,TX},\n    tsteps::Union{StepRangeLen,AbstractVector},\n    mu,\n    iv::Type{AbstractIndependantVariable};\n    solver = Vern9(),\n    reltol = 1e-14,\n    abstol = 1e-14,\n) where TX\n\nPropagate the state x0 from tsteps[1] to tsteps[end], returning the solution at each time in tsteps. The final argument can be specified to set the desired independant variable (i.e., Time or ArcLength).\n\nArguments\n\nx0::SVector{6,TX}: The initial state with x0 = [r0; v0].\ntsteps::Union{StepRangeLen,AbstractVector}: The time steps at which a state x(t) is desired.\nmu::Real: The mass parameter of the three-body system (mu = m2 / (m1 + m2)).\niv::Type{AbstractIndependantVariable}: The desired independant variable.\n\nKeywords\n\nsolver::OrdinaryDiffEq.AbstractODESolver: The solver to use for the integration.\nreltol: The relative tolerance for the solver.\nabstol: The absolute tolerance for the solver.\n\nReturns\n\nVector{SVector{6,TX}}: A vector of states at each time in tsteps.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#CRTBPNaturalMotion.propagate_return_arc_length-Union{Tuple{TT}, Tuple{TX}, Tuple{SVector{6, TX}, Tuple{TT, TT}, Any}} where {TX, TT}","page":"Public API","title":"CRTBPNaturalMotion.propagate_return_arc_length","text":"propagate_return_arc_length(\n    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu;\n    solver = Vern9(),\n    reltol = 1e-14,\n    abstol = 1e-14,\n) where {TX,TT}\n\nPropagate the state x0 over the time span tspan and return the arclength along the trajectory. The final argument can be specified to set the desired independant variable (i.e., Time or ArcLength).\n\narguments\n\nx0::SVector{6,TX}: the initial state with x0 = [r0; v0].\ntspan::Tuple{TT,TT}: the time span forwhich to solve the ode.\nmu::Real: the mass parameter of the three-body system (mu = m2 / (m1 + m2)).\n\nkeywords\n\nsolver::OrdinaryDiffEq.AbstractODESolver: The solver to use for the integration.\nreltol: The relative tolerance for the solver.\nabstol: The absolute tolerance for the solver.\n\nreturns\n\nTX: The arc-length along the curve\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#CRTBPNaturalMotion.propagate_return_final_state-Union{Tuple{TT}, Tuple{TX}, Tuple{SVector{6, TX}, Tuple{TT, TT}, Any, Function, Type{Time}}} where {TX, TT}","page":"Public API","title":"CRTBPNaturalMotion.propagate_return_final_state","text":"propagate_return_final_state(\n    x0::SVector{6,TX},\n    tspan::Tuple{TT,TT},\n    mu,\n    term_cond::Function [, iv::Type{AbstractIndependantVariable}];\n    solver = Vern9(),\n    reltol = 1e-14,\n    abstol = 1e-14,\n) where {TX,TT}\n\nPropagate the state x0 over the time span tspan and return the final state. At the first point at which term_cond is zero (if any), the integration will terminate. The final argument can be specified to set the desired independant variable (i.e., Time or ArcLength).\n\nArguments\n\nx0::SVector{6,TX}: The initial state with x0 = [r0; v0].\ntspan::Tuple{TT,TT}: The span of the independant variable forwhich to solve the ODE.\nmu::Real: The mass parameter of the three-body system (mu = m2 / (m1 + m2)).\nterm_cond::Function: The DifferentialEquations.jl continuous callback condition function   that will terminate the integration when it is zero. Function should be of the form   term_cond(x,t,integ), where x is the state, t is the independant variable, and   integ is the integrator struct (see DifferentialEquations.jl callback documentation).\niv::Type{AbstractIndependantVariable}: The desired independant variable (Time if   unspecified).\n\nKeywords\n\nsolver::OrdinaryDiffEq.AbstractODESolver: The solver to use for the integration.\nreltol: The relative tolerance for the solver.\nabstol: The absolute tolerance for the solver.\n\nReturns\n\nSVector{6,TX}: The final state of the ODE solution.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#CRTBPNaturalMotion.propagate_return_final_state-Union{Tuple{TT}, Tuple{TX}, Tuple{SVector{6, TX}, Tuple{TT, TT}, Any, Type{Time}}} where {TX, TT}","page":"Public API","title":"CRTBPNaturalMotion.propagate_return_final_state","text":"propagate_return_final_state(\n    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu [, iv::Type{AbstractIndependantVariable}];\n    solver = Vern9(),\n    reltol = 1e-14,\n    abstol = 1e-14,\n) where {TX,TT}\n\nPropagate the state x0 over the time span tspan and return the final state. The final argument can be specified to set the desired independant variable (i.e., Time or ArcLength).\n\narguments\n\nx0::SVector{6,TX}: The initial state with x0 = [r0; v0].\ntspan::Tuple{TT,TT}: The span of the independant variable forwhich to solve the ode.\nmu::Real: The mass parameter of the three-body system (mu = m2 / (m1 + m2)).\niv::Type{AbstractIndependantVariable}: The desired independant variable (Time if   unspecified).\n\nkeywords\n\nsolver::OrdinaryDiffEq.AbstractODESolver: The solver to use for the integration.\nreltol: The relative tolerance for the solver.\nabstol: The absolute tolerance for the solver.\n\nreturns\n\nSVector{6,TX}: The final state of the ode solution.\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#CRTBPNaturalMotion.propagate_return_final_stm-Union{Tuple{TT}, Tuple{TX}, Tuple{SVector{6, TX}, Tuple{TT, TT}, Any}} where {TX, TT}","page":"Public API","title":"CRTBPNaturalMotion.propagate_return_final_stm","text":"propagate_return_final_stm(\n    x0::SVector{6,TX}, tspan::Tuple{TT,TT}, mu;\n    solver = Vern9(),\n    reltol = 1e-14,\n    abstol = 1e-14,\n) where {TX,TT}\n\nPropagate the state x0 over the time span tspan and return the final state transition matrix (STM).\n\narguments\n\nx0::SVector{6,TX}: The initial state with x0 = [r0; v0].\ntspan::Tuple{TT,TT}: The span of the independant variable forwhich to solve the ode.\nmu::Real: The mass parameter of the three-body system (mu = m2 / (m1 + m2)).\n\nkeywords\n\nsolver::OrdinaryDiffEq.AbstractODESolver: The solver to use for the integration.\nreltol: The relative tolerance for the solver.\nabstol: The absolute tolerance for the solver.\n\nreturns\n\nSMatrix{6,6,TX,36}: The stm for the full trajectory\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#CRTBPNaturalMotion.propagate_return_time-Union{Tuple{TT}, Tuple{TX}, Tuple{SVector{6, TX}, Tuple{TT, TT}, Any}} where {TX, TT}","page":"Public API","title":"CRTBPNaturalMotion.propagate_return_time","text":"propagate_return_time(\n    x0::SVector{6,TX}, sspan::Tuple{TT,TT}, mu;\n    solver = Vern9(),\n    reltol = 1e-14,\n    abstol = 1e-14,\n) where {TX,TT}\n\nPropagate the state x0 over the arc-length span sspan and return the arclength along the trajectory. The final argument can be specified to set the desired independant variable (i.e., Time or ArcLength).\n\narguments\n\nx0::SVector{6,TX}: the initial state with x0 = [r0; v0].\nsspan::Tuple{TT,TT}: the arc-length span forwhich to solve the ode.\nmu::Real: the mass parameter of the three-body system (mu = m2 / (m1 + m2)).\n\nkeywords\n\nsolver::OrdinaryDiffEq.AbstractODESolver: The solver to use for the integration.\nreltol: The relative tolerance for the solver.\nabstol: The absolute tolerance for the solver.\n\nreturns\n\nTX: The time along the curve\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#Utility-Functions","page":"Public API","title":"Utility Functions","text":"","category":"section"},{"location":"lib/public/","page":"Public API","title":"Public API","text":"Modules = [CRTBPNaturalMotion]\nPages = [\"utils.jl\"]\nPrivate = false","category":"page"},{"location":"lib/public/#CRTBPNaturalMotion.jacobi_integral-Tuple{Any, Any}","page":"Public API","title":"CRTBPNaturalMotion.jacobi_integral","text":"jacobi_integral(x, μ)\n\nComputes the Jacobi integral for the CRTBP given the state x and mass parameter μ.\n\nArguments\n\nx::AbstractVector: State vector.\nμ::Real: Mass parameter.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = CRTBPNaturalMotion","category":"page"},{"location":"#CRTBPNaturalMotion","page":"Home","title":"CRTBPNaturalMotion","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Circular Restricted Three-Body Problem (CRTBP)NaturalMotion is a package for computing periodic orbits and they're associated stable/unstable invariant manifolds, as well as performing various analyses in the CRTBP model assuming only natural motion (i.e., no control)","category":"page"},{"location":"#Simple-Example-of-Computing-Halo-Orbit","page":"Home","title":"Simple Example of Computing Halo Orbit","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We'll first define some initial guess for a reference state of a periodic orbit of Type A. For this class of periodic orbits, we'll need to specify a position component in the hatx- and hatz-directions of the synodic (rotating) reference frame (i.e., r_x and r_z), as well as a velocity component in the haty-direction (v_y). This will correspond to the full state mathbfx = r_x 00 r_z 00 v_y 00^T, which lies in the hatx-hatz plane with velocity perpendicular to this plane. Note we'll also define the mass parameter for the Earth-Moon system, which is given by mu = m_2  (m_1 + m_2) (where here, m_1 is the mass of the Earth and m_2 is the mass of the Moon).","category":"page"},{"location":"","page":"Home","title":"Home","text":"# CRTBP mass parameter\nmu = 1.21506038e-2\n\n# Define initial state guess\nrx_guess = 0.8233851820\nrz_guess = -0.0222775563\nvy_guess = 0.1341841703\nnothing # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"We'll now compute the Jacobi integral (the only constant of motion for the CRTBP model) that corresponds to this initial state. Note we'll also load CRTBPNaturalMotion.jl and the StaticArrays.jl package so we can use their fast arrays here (and later).","category":"page"},{"location":"","page":"Home","title":"Home","text":"using CRTBPNaturalMotion\nusing StaticArrays\n\nC_guess = jacobi_integral(SA[rx_guess, 0.0, rz_guess, 0.0, vy_guess, 0.0], mu)\nnothing # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, we can now construct a TypeAPeriodicOrbit employing our initial guess, which will get corrected to a periodic orbit that satisfies an extra constraint (as a keyword argument). We'll chose a constraint that requires the Jacobi integral of the final orbit to be close to, but not quite the same, as the Jacobi integral corresponding to our initial guess. Note we'll also specify the keyword argument N_cross, which sets the number of times the periodic orbit will pass through the hatx-hatz plane in half and orbital period (i.e., N_cross = 1 means we'll pass through this plane twice if traveling around the full orbit).","category":"page"},{"location":"","page":"Home","title":"Home","text":"halo = TypeAPeriodicOrbit(\n    rx_guess, rz_guess, vy_guess, mu;\n    N_cross = 1,\n    constraint = (:jacobi_integral, C_guess - 0.001),\n)","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, we have some useful functions for using this TypeAPeriodicOrbit. For now, we'll just use the get_full_orbit function to obtain a DifferentialEquations.jl solution containing the states around the full orbit so we can plot using CairoMakie.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using CairoMakie\n\ntraj = get_full_orbit(halo)\n\nfig = Figure(); \nax  = Axis3(\n    fig[1,1]; \n    aspect = :data,\n    xlabel = L\"$r_x$, LU\",\n    ylabel = L\"$r_y$, LU\",\n    zlabel = L\"$r_z$, LU\",\n)\nlines!(ax, traj[1,:], traj[2,:], traj[3,:])\nfig # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"Let's now try to compute a bunch of periodic orbits using a natural continuation based approach, where we'll gradually adjust the value of the constraint while solving for new periodic orbits repeatedly (taking the previous solution as the next initial guess). Note that there are much better methods of performing continuation to compute families of periodic orbit, but this is quick example that is simple to implement.","category":"page"},{"location":"","page":"Home","title":"Home","text":"# Define starting orbit\nstart_orbit = TypeAPeriodicOrbit(\n    rx_guess, rz_guess, vy_guess, mu;\n    N_cross = 1,\n    constraint = (:jacobi_integral, C_guess),\n)\n\n# Set desired rx and number of steps to take\nC_final = C_guess - 0.05\nn_steps = 100\n\n# Perform continuation in loop\nCs = range(C_guess, C_final; length = n_steps)\norbits = [start_orbit]\nfor (i, C) in enumerate(Cs)\n    new_orbit = TypeAPeriodicOrbit(\n        orbits[end];\n        constraint = (:jacobi_integral, C)\n    )\n    push!(orbits, new_orbit)\n\n    # Update our figure\n    if mod(i, 10) == 0\n        new_traj = get_full_orbit(new_orbit)\n        lines!(ax, new_traj[1,:], new_traj[2,:], new_traj[3,:]; color = :red)\n    end\nend\nfig # hide","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
