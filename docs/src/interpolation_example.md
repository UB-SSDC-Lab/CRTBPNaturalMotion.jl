
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

Next, we'll create two different `FastChebInterpolation`s, one representing an interpolation of the periodic orbit, and another representing an approximation of the orbit employing a least-squares fit to a Chebyshev polynomial of the same order. Note that we'll generate both as functions of the scaled `ArcLength` about one full period of the orbit, i.e., `τ` ``\in [0, 1]``.

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

Clearly, we get similar values from either polynomial representation. Note we can also easily compute the derivative of the polynomial approximations with: 
```@example interpolation
    derivative(orbit_lsqf, τ) 
```
and second derivative with:
```@example interpolation
    second_derivative(orbit_lsqf, τ)
```

Let's now see the approximation error using both methods.

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
All `InvariantManifold`s in `CRTBPNaturalMotion.jl` can also be approximated with Chebyshev polynomials. We can construct a bi-variate approximation of a full manifold surface or a uni-variate approximation of a manifold cross-section (both with interpolation and least-squares regression). To demonstrate this, we'll first plot trajectories that lie on the stable invariant manifold of an ``L_1`` halo orbit as follows:

```@example interpolation
    # Define halo orbit and manifold
    halo = get_jpl_orbits(;
        sys         = "earth-moon",
        family      = "halo",
        libr        = 1, 
        jacobimin   = 3.104,
        jacobimax   = 3.105,
    )[1]
    man = InvariantManifold(halo, 50.0 / halo.DU)

    # Some convenience variables for computing manifold
    man_start = 0.75 
    man_length = 6e5 / halo.DU

    # Lets plot some trajectories on the manifold
    fig = Figure()
    ax  = Axis(fig[1,1]; aspect = DataAspect(), xlabel = L"$r_x$, DU", ylabel = L"$r_y$, DU")
    τs  = range(0.0, 1.0; length=16)[1:end-1]
    for τ in τs
        # Get the initial manifold trajectory for better visualization of what is going on
        init_man_traj = get_stable_manifold_trajectory(
            man, true, τ, ArcLength, 0.0, 1.0, ArcLength;
            manifold_length = man_start 
        )

        # Get the manifold trajectory that we'll actually use
        man_traj = get_stable_manifold_trajectory(
            man, true, τ, ArcLength, man_start, 1.0, ArcLength;
            manifold_length = man_length # 8e5 km divided by the distance unit
        )
        lines!(ax, init_man_traj[1,:], init_man_traj[2,:]; color = (:black, 0.5))
        lines!(ax, man_traj[1,:], man_traj[2,:]; color = :blue)
    end
    fig # hide
```

Note that the blue portions of each trajectory correspond to the part of the stable invariant manifold that will be considered while the grey portions will be ignored, i.e., we'll ignore part of the manifold that lies very close to the halo orbit (i.e., the last 0.75 DU in terms of the `ArcLength` along the manifold), and will consider a portion of the manifold with a length of ``6 \times 10^5`` km of the stable manifold in terms of the `ArcLength`. As such, `τ2 = 1.0` corresponds to states on the manifold that are a distance of ``6 \times 10^5 + 0.75 \times`` `halo.DU` km, in terms of `ArcLength` along the manifold. 

### Approximation of a Manifold Cross-Section

Let's first consider approximating a single cross section of this stable manifold. We'll consider a cross section that is located ``5 \times 10^5`` km from the halo orbit, in terms of the `ArcLength` along the manifold. To visualize this cross-section, we can plot it as follows:

```@example interpolation
    fig = Figure()
    ax = Axis(fig[1,1]; aspect = DataAspect(), xlabel = L"$r_x$, DU", ylabel = L"$r_y$, DU")
    
    # First plot the halo orbit itself
    halo_traj = get_full_orbit(halo)
    lines!(ax, halo_traj[1,:], halo_traj[2,:], halo_traj[3,:]; color = :black)

    # Define cross-section location
    cross_section_location = 3e5 / halo.DU
    cross_section_τ2 = cross_section_location / man_length

    # Plot cross-section
    N = 100
    τ1s = range(0.0, 1.0; length = N)
    cs_states = zeros(6, N)
    for (i, τ1) in enumerate(τ1s)
        manifold_state = get_stable_manifold_trajectory(
            man, true, τ1, ArcLength, man_start, cross_section_τ2, ArcLength;
            manifold_length = man_length, 
        )[end]
        cs_states[:,i] .= manifold_state
    end
    lines!(ax, cs_states[1,:], cs_states[2,:], cs_states[3,:]; color = :blue)

    # Plot manifold trajectories
    τs  = range(0.0, 1.0; length=16)[1:end-1]
    for τ in τs
        man_traj = get_stable_manifold_trajectory(
            man, true, τ, ArcLength, 0.0, 1.0, ArcLength;
            manifold_length = man_start + man_length
        )
        lines!(ax, man_traj[1,:], man_traj[2,:]; color = (:black, 0.5))
    end

    fig # hide
```

The Chebyshev polynomial approximations can be constructed with:
```@example interpolation
    # First define polynomial order and the number of points to use in the 
    # least-squares fit
    τ1_order = 120
    τ1_npoints = 200


    man_5e5_cs_interp = generate_stable_manifold_cross_section_cheb_interpolant(
        man, true, τ1_order, ArcLength, cross_section_τ2, ArcLength;
        start_cond = man_start,
        manifold_length = man_length,
    )
    man_5e5_cs_lsqf = generate_stable_manifold_cross_section_cheb_approximation(
        man, true, τ1_npoints, τ1_order, ArcLength, cross_section_τ2, ArcLength;
        start_cond = man_start,
        manifold_length = man_length,
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
            man, true, τ1, ArcLength, man_start, cross_section_τ2, ArcLength;
            manifold_length = man_length, 
        )[end]
        interp_error[:,i] .= value(man_5e5_cs_interp, τ1) - manifold_state
        lsqf_error[:,i] .= value(man_5e5_cs_lsqf, τ1) - manifold_state
    end

    # Plot position and velocity errors in x-direction with units of m and m/s
    fig = Figure()
    ax_rx = Axis(fig[1,1]; xlabel = L"\tau_1", ylabel=L"$r_x$ error (m)")
    lines!(
        ax_rx, τ1_check_steps, lsqf_error[1,:]*(halo.DU*1000); 
        linewidth   = 1,
        color       = :green, 
        label       = "Least-Squares approx.",
    )
    lines!(
        ax_rx, τ1_check_steps, interp_error[1,:]*(halo.DU*1000); 
        linewidth   = 1,
        color       = :blue,
        label       = "Interpolation",
    )
    axislegend(ax_rx)

    ax_vx = Axis(fig[2,1]; xlabel = L"\tau_1", ylabel=L"$v_x$ error (m/s)")
    lines!(
        ax_vx, τ1_check_steps, lsqf_error[4,:]*(halo.DU*1000/halo.TU); 
        linewidth   = 1,
        color       = :green, 
        label       = "Least-Squares approx.",
    )
    lines!(
        ax_vx, τ1_check_steps, interp_error[4,:]*(halo.DU*1000/halo.TU); 
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
    τ1_order = 80
    τ2_order = 80
    τ1_npionts = 150
    τ2_npoints = 150

    # Construct the interpolant
    man_interp = generate_stable_manifold_cheb_interpolant(
        man, true, τ1_order, ArcLength, τ2_order, ArcLength;
        start_cond = man_start,
        manifold_length = man_length, 
    )

    # Compute manifold interpolation error
    N = 1000
    τ1_check_steps = range(0.0, 1.0, length = N)
    τ2_check_steps = range(0.0, 1.0, length = N)
    interp_errors  = zeros(6, N, N)
    for (i,τ1) in enumerate(τ1_check_steps)
        manifold_states = get_stable_manifold_trajectory(
            man, true, τ1, ArcLength, man_start, τ2_check_steps, ArcLength;
            manifold_length = man_length,
        )
        for j in eachindex(manifold_states)
            τ2 = τ2_check_steps[j]
            prop_val = manifold_states[j]
            interp_val = value(man_interp, τ1, τ2)
            interp_errors[:,i,j] .= interp_val - prop_val
        end
    end

    # Plot interpolation error surface for x position and velocity
    fig1 = Figure()
    ax1 = Axis3(fig1[1,1]; xlabel=L"\tau_1", ylabel=L"\tau_2", zlabel=L"$r_x$ error, m")
    surface!(
        ax1, τ1_check_steps, τ2_check_steps, interp_errors[1,:,:]*(halo.DU*1000.0); 
        colormap = :viridis,
    )
    fig1 # hide
```
