
# Interpolation in CRTBPNaturalMotion.jl
`CRTBPNaturalMotion.jl` implements several methods for approximating periodic orbits or their associated stable/unstable invariant manifolds. Currently all methods employ Chebyshev polynomials from [`FastChebInterp.jl`](https://github.com/JuliaMath/FastChebInterp.jl).

## Interpolation of Periodic Orbits
All periodic orbits in `CRTBPNaturalMotion.jl` can be approximated via Chebyshev polynomials, either through interpolation or least-squares regression. To demonstrate this functionality, we'll employ an axial periodic orbit about ``L_4`` in the Earth-Moon system:

```@example interpolation
    using CRTBPNaturalMotion

    axial_orbit = get_jpl_orbits(;
        sys         = "earth-moon",
        family      = "axial",
        libr        = 4,
        periodmin   = 20.0,
        periodmax   = 21.0,
        periodunits = "d",
    )[1]
```

Next, we'll create two different `FastChebInterpolation`s, one representing an interpolation of the periodic orbit, and another representing an approximation of the orbit employing a least-squares fit to a Chebyshev polynomial of the same order. Note that we'll generate both as functions of the scaled `ArcLength`, i.e., `τ` ``\in [0, 1]``.

```@example interpolation
    # Define order of both polynomials
    polynomial_order = 100

    # Define number of points to employ in least-squares fit
    lsq_npoints = 150

    # Generate interpolation of orbit
    orbit_interp = generate_periodic_orbit_cheb_interpolation(
        axial_orbit, polynomial_order, ArcLength,
    )

    # Generate approximation with least-squares regression
    orbit_lsqf = generate_periodic_orbit_cheb_approximation(
        axial_orbit, lsq_npoints, polynomial_order, ArcLength,
    )

    # Let's get the value of each polynomial for some τ
    τ = 0.3256
    println(value(orbit_interp, τ))
    println(value(orbit_lsqf, τ))
```

Clearly, we get similar values from either polynomial representation. Let's now see the approximation error using both methods.

```@example interpolation
    using CairoMakie

    N = 10000
    τs = range(0.0, 1.0; length = N)
    interp_errors = zeros(6, N)
    lsqf_errors = zeros(6, N)
    for (i, τ) in enumerate(τs)
        # Propagate dynamics from initial reference state (uses Vern9() 
        # with absolute and relative tolerance of 1e-14 be default)
        prop_state = CRTBPNaturalMotion.propagate_from_initial_conditions(
            axial_orbit, τ, ArcLength,
        )

        # Compute errors
        interp_errors[:, i] .= prop_state - value(orbit_interp, τ)
        lsqf_errors[:, i] .= prop_state - value(orbit_lsqf, τ)
    end

    # Plot position and velocity errors in x-direction with units of m
    fig = Figure()
    ax_rx = Axis(fig[1,1]; xlabel = L"\tau", ylabel=L"$r_x$ error (m)")
    lines!(
        ax_rx, τs, lsqf_errors[1,:]*(axial_orbit.DU*1000); 
        linewidth   = 1,
        color       = :green, 
        label       = "Least-Squares approx.",
    )
    lines!(
        ax_rx, τs, interp_errors[1,:]*(axial_orbit.DU*1000); 
        linewidth   = 1,
        color       = :blue,
        label       = "Interpolation",
    )
    axislegend(ax_rx)

    ax_vx = Axis(fig[2,1]; xlabel = L"\tau", ylabel=L"$v_x$ error (m/s)")
    lines!(
        ax_vx, τs, lsqf_errors[4,:]*(axial_orbit.DU*1000/axial_orbit.TU); 
        linewidth   = 1,
        color       = :green, 
        label       = "Least-Squares approx.",
    )
    lines!(
        ax_vx, τs, interp_errors[4,:]*(axial_orbit.DU*1000/axial_orbit.TU); 
        linewidth   = 1,
        color       = :blue,
        label       = "Interpolation",
    )
    fig # hide
```

As can be seen, both approximations are quite good, with the least-squares approximation exhibiting the smallest maximum error whereas the interpolation exhibits smaller error for values of the parameter ``\tau`` near the boundaries.

## Stable Invariant Manifold
All `InvariantManifold`s in `CRTBPNaturalMotion.jl` can also be approximated with Chebyshev polynomials. We can construct a bi-variate approximation of a full manifold surface or a uni-variate approximation of a manifold cross-section (both with interpolation and least-squares regression). To demonstrate this, we'll first construct an `InvariantManifold` for an ``L_1`` halo orbit as follows:

```@example interpolation
    # Define units for convenience
    G       = 6.673e-20
    m_1     = 5.9742e24         # kg
    m_2     = 7.3483e22         # kg
    r_12    = 384400.0          # km
    DU      = 384400.0
    TU      = 1.0 / sqrt((G*(m_1 + m_2)) / DU^3)

    # Define halo orbit and manifold
    halo = TypeAPeriodicOrbit(
        0.8233851820, -0.0222775563, 0.1341841703, 1.21506038e-2, TU, DU;
        N_cross    = 1,
        constraint = (:x_start_coordinate, 0.8233851820),
        ode_solver = Vern9(),
        ode_reltol = 1e-14,
        ode_abstol = 1e-14,
    )
    man = InvariantManifold(halo, 1e-3 / 384400.0)

    # Lets plot some trajectories on the manifold
    fig = Figure()
    ax  = Axis(fig[1,1]; aspect = DataAspect())
    τs  = range(0.0, 1.0; length=16)[1:end-1]
    for τ in τs
        man_traj = get_stable_manifold_trajectory(
            man, true, τ, ArcLength, 0.0, 1.0, ArcLength;
            manifold_length = 8e5 / DU # 8e5 km divided by the distance unit
        )
        lines!(ax, man_traj[1,:], man_traj[2,:]; color = :blue)
    end
    fig # hide
```
Note that we're specifying that we're only interested in the last ``8 \times 10^5`` km of the stable manifold, such that `τ2 = 1.0` corresponds states on the manifold that are a distance of ``8 \times 10^5`` km, in terms of `ArcLength` along the manifold. 

### Approximation of a Manifold Cross-Section

Let's first consider approximating a single cross section of this stable manifold. We'll consider a cross section that is located ``5 \times 10^5`` km from the halo orbit, in terms of the `ArcLength` along the manifold. To visualize this cross-section, we can plot it as follows:

```@example interpolation
    fig = Figure()
    ax = Axis(fig[1,1]; aspect = DataAspect())
    
    # First plot the halo orbit itself
    halo_traj = get_full_orbit(halo)
    lines!(ax, halo_traj[1,:], halo_traj[2,:], halo_traj[3,:]; color = :black)

    # Plot cross-section
    N = 100
    τ1s = range(0.0, 1.0; length = N)
    cs_states = zeros(6, N)
    for (i, τ1) in enumerate(τ1s)
        manifold_state = get_stable_manifold_trajectory(
            man, true, τ1, ArcLength, 0.0, 5e5 / 8e5, ArcLength;
            manifold_length = 8e5 / DU, 
        )[end]
        cs_states[:,i] .= manifold_state
    end
    lines!(ax, cs_states[1,:], cs_states[2,:], cs_states[3,:]; color = :blue)

    fig # hide
```


The Chebyshev polynomial approximations can be constructed with:
```@example interpolation
    # First define polynomial order and the number of points to use in the 
    # least-squares fit
    τ1_order = 120
    τ1_npoints = 200


    man_5e5_cs_interp = generate_stable_manifold_cross_section_cheb_interpolant(
        man, true, τ1_order, ArcLength, (5e5 / 8e5), ArcLength;
        start_cond = 0.0,
        manifold_length = 8e5 / DU,
    )
    man_5e5_cs_lsqf = generate_stable_manifold_cross_section_cheb_approximation(
        man, true, τ1_npoints, τ1_order, ArcLength, 5e5 / 8e5, ArcLength;
        start_cond = 0.0,
        manifold_length = 8e5 / DU,
    )

    τ1 = 0.3256
    println(value(man_5e5_cs_interp, τ1))
    println(value(man_5e5_cs_lsqf, τ1))
```

Again, lets look at the approximation error using both methods by plotting with `CairoMakie.jl`:
```@example interpolation
    N = 1000
    τ1_check_steps = range(0.0, 1.0; length = N)
    interp_error = zeros(6, N)
    lsqf_error = zeros(6, N)
    for (i, τ1) in enumerate(τ1_check_steps)
        manifold_state = get_stable_manifold_trajectory(
            man, true, τ1, ArcLength, 0.0, 5e5 / 8e5, ArcLength;
            manifold_length = 8e5 / 384400.0, 
        )[end]
        interp_error[:,i] .= value(man_5e5_cs_interp, τ1) - manifold_state
        lsqf_error[:,i] .= value(man_5e5_cs_lsqf, τ1) - manifold_state
    end

    # Plot position and velocity errors in x-direction with units of m and m/s
    fig = Figure()
    ax_rx = Axis(fig[1,1]; xlabel = L"\tau_1", ylabel=L"$r_x$ error (m)")
    lines!(
        ax_rx, τ1_check_steps, lsqf_error[1,:]*(DU*1000); 
        linewidth   = 1,
        color       = :green, 
        label       = "Least-Squares approx.",
    )
    lines!(
        ax_rx, τ1_check_steps, interp_error[1,:]*(DU*1000); 
        linewidth   = 1,
        color       = :blue,
        label       = "Interpolation",
    )
    axislegend(ax_rx)

    ax_vx = Axis(fig[2,1]; xlabel = L"\tau_1", ylabel=L"$v_x$ error (m/s)")
    lines!(
        ax_vx, τ1_check_steps, lsqf_error[4,:]*(DU*1000/TU); 
        linewidth   = 1,
        color       = :green, 
        label       = "Least-Squares approx.",
    )
    lines!(
        ax_vx, τ1_check_steps, interp_error[4,:]*(DU*1000/TU); 
        linewidth   = 1,
        color       = :blue,
        label       = "Interpolation",
    )
    fig # hide
```

### Approximation of a Full Manifold Surface
Finally, we'll compute an approximation of the full `InvariantManifold` surface. In the following code example, we'll construct an interpolant of the full manifold and will compute the approximation error.

```@example interpolation
    # First define order of polynomial for each variable and the number of points
    # in each dimension to use in the least-squares fit
    τ1_order = 150
    τ2_order = 150
    τ1_npionts = 400
    τ2_npoints = 400

    # Construct the interpolant
    man_interp = generate_stable_manifold_cheb_interpolant(
        man, true, τ1_order, ArcLength, τ2_order, ArcLength;
        start_cond = 0.0 / DU,
        manifold_length = 8e5 / DU, 
    )

    # This is how we'd construct the least-squares fit, but it takes a long 
    # time to run so we won't run it now...
    # man_lsqf = generate_stable_manifold_cheb_approximation(
    #     man, true, τ1_npoints, τ1_order, ArcLength, τ2_npoints, τ2_order, ArcLength;
    #     start_cond = 0.0,
    #     manifold_length = 8e5 / DU, 
    # )

    # Compute manifold interpolation error
    N = 1000
    τ1_check_steps = range(0.0, 1.0, length = N)
    τ2_check_steps = range(0.0, 1.0, length = N)
    errors         = zeros(6, N, N)
     for (i,τ1) in enumerate(τ1_check_steps)
        manifold_states = get_stable_manifold_trajectory(
            man, true, τ1, ArcLength, 0.0 / DU, τ2_check_steps, ArcLength;
            manifold_length = 8e5 / DU,
        )
        for j in eachindex(manifold_states)
            τ2 = τ2_check_steps[j]
            prop_val = manifold_states[j]
            interp_val = value(man_interp, τ1, τ2)
            errors[:,i,j] .= interp_val - prop_val
        end
    end

    # Plot interpolation error surface for x position and velocity
    fig = Figure()
    ax1 = Axis3(fig[1,1]; xlabel=L"\tau_1", ylabel=L"\tau_2", zlabel=L"$r_x$ error, m")
    surface!(
        ax1, τ1_check_steps, τ2_check_steps, errors[1,:,:]*(DU*1000.0); 
        colormap = :viridis,
    )
    fig #hide
```
