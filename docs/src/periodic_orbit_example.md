# Computing Halo Orbits
We'll first define some initial guess for a reference state of a periodic orbit of Type A. For this class of periodic orbits, we'll need to specify a position component in the ``\hat{x}``- and ``\hat{z}``-directions of the synodic (rotating) reference frame (i.e., ``r_x`` and ``r_z``), as well as a velocity component in the ``\hat{y}``-direction (``v_y``). This will correspond to the full state ``\mathbf{x} = [r_x, 0.0, r_z, 0.0, v_y, 0.0]^T``, which lies in the ``\hat{x}-\hat{z}`` plane with velocity perpendicular to this plane. Note we'll also define the mass parameter for the Earth-Moon system, which is given by ``\mu = m_2 / (m_1 + m_2)`` (where here, ``m_1`` is the mass of the Earth and ``m_2`` is the mass of the Moon).

```@example simple_halo
# CRTBP mass parameter
mu = 1.21506038e-2

# Define initial state guess
rx_guess = 0.8233851820
rz_guess = -0.0222775563
vy_guess = 0.1341841703
nothing # hide
```

We'll now compute the Jacobi integral (the only constant of motion for the CRTBP model) that corresponds to this initial state. Note we'll also load CRTBPNaturalMotion.jl and the StaticArrays.jl package so we can use their fast arrays here (and later).
```@example simple_halo
using CRTBPNaturalMotion
using StaticArrays

C_guess = jacobi_integral(SA[rx_guess, 0.0, rz_guess, 0.0, vy_guess, 0.0], mu)
nothing # hide
```

Finally, we can now construct a `TypeAPeriodicOrbit` employing our initial guess, which will get corrected to a periodic orbit that satisfies an extra `constraint` (as a keyword argument). We'll chose a constraint that requires the Jacobi integral of the final orbit to be close to, but not quite the same, as the Jacobi integral corresponding to our initial guess. Note we'll also specify the keyword argument `N_cross`, which sets the number of times the periodic orbit will pass through the ``\hat{x}-\hat{z}`` plane in half and orbital period (i.e., `N_cross = 1` means we'll pass through this plane twice if traveling around the full orbit).
```@example simple_halo
halo = TypeAPeriodicOrbit(
    rx_guess, rz_guess, vy_guess, mu;
    N_cross = 1,
    constraint = (:jacobi_integral, C_guess - 0.001),
)
```

Finally, we have some useful functions for using this `TypeAPeriodicOrbit`. For now, we'll just use the `get_full_orbit` function to obtain a DifferentialEquations.jl solution containing the states around the full orbit so we can plot using CairoMakie.jl.
```@example simple_halo
using CairoMakie

traj = get_full_orbit(halo)

fig = Figure(); 
ax  = Axis3(
    fig[1,1]; 
    aspect = :data,
    xlabel = L"$r_x$, DU",
    ylabel = L"$r_y$, DU",
    zlabel = L"$r_z$, DU",
)
lines!(ax, traj[1,:], traj[2,:], traj[3,:])
fig # hide
```

Let's now try to compute a bunch of periodic orbits using a natural continuation based approach, where we'll gradually adjust the value of the `constraint` while solving for new periodic orbits repeatedly (taking the previous solution as the next initial guess). Note that there are much better methods of performing continuation to compute *families* of periodic orbit, but this is quick example that is simple to implement.
```@example simple_halo
# Define starting orbit
start_orbit = TypeAPeriodicOrbit(
    rx_guess, rz_guess, vy_guess, mu;
    N_cross = 1,
    constraint = (:jacobi_integral, C_guess),
)

# Set desired rx and number of steps to take
C_final = C_guess - 0.05
n_steps = 100

# Perform continuation in loop
Cs = range(C_guess, C_final; length = n_steps)
orbits = [start_orbit]
for (i, C) in enumerate(Cs)
    new_orbit = TypeAPeriodicOrbit(
        orbits[end];
        constraint = (:jacobi_integral, C)
    )
    push!(orbits, new_orbit)

    # Update our figure
    if mod(i, 10) == 0
        new_traj = get_full_orbit(new_orbit)
        lines!(ax, new_traj[1,:], new_traj[2,:], new_traj[3,:]; color = :red)
    end
end
fig # hide
```
