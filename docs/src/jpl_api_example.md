# Querying JPL's CRTBP Poincare Catalog of Periodic Orbits
`CRTBPNaturalMotion.jl` provides an interface to JPL's CRTBP Poincare catalog of periodic orbits through the exported function `get_jpl_orbits`. Through several keyword arguments (see the [`get_jpl_orbits`](@ref) docstring), which directly correspond to the in JPL's [web application](https://ssd.jpl.nasa.gov/tools/periodic_orbits.html) and [API documentation](https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html), an `OrbitSet` of periodic orbits can be obtained for various three-body systems (e.g., Earth-Moon, Sun-Earth, Saturn-Enceladus, etc...).

## Earth-Moon Halo orbits

As a simple example, we can obtain the full set of periodic orbits about the ``L_1`` Libration point as follows:
```@example jpl_api
using CRTBPNaturalMotion

orbits = get_jpl_orbits(;
    sys     = "earth-moon",
    family  = "halo",
    libr    = 1,
)
```
As can be seen, a total of 5731 orbits are contained within the `OrbitSet` exhibiting a range of different orbital periods, Jacobi constants, and stability values. 

Perhaps we're only interested in halo orbits exhibiting a stability index that is less than 5, with orbital periods between 10 and 11 days. We can refine our query using some additional keyword arguments as follows:
```@example jpl_api
stable_orbits = get_jpl_orbits(;
    sys         = "earth-moon",
    family      = "halo",
    libr        = 1,
    periodmin   = 10.0,
    periodmax   = 11.0,
    periodunits = "d",
    stabmax     = 5.0,
)
```
The returned `OrbitSet` is implemented as a sub-type of Julia's `AbstractVector{T}` and therefore functions such as `length(::OrbitSet)` and indexing, i.e., 
```@example jpl_api
first_orbit = stable_orbits[1]
```
and 
```@example jpl_api
stable_orbits_subset = stable_orbits[[3,5,10,12]]
```
can be employed to obtain a single `GeneralPeriodicOrbit` or a subset of the original `OrbitSet`. Note that an `OrbitSet` is essentially a `Vector{GeneralPeriodicOrbit}` with some extra information about the three-body `System` and the contained periodic orbits. 

We can obtain a full trajectory for the periodic orbit as shown in [Computing Halo Orbits](@ref) using the `get_full_orbit` method. For example, several orbits within `stable_orbits` can be plotted using `CairoMakie.jl` as follows:
```@example jpl_api
using CairoMakie

fig = Figure(; size = (700, 1000))
ax  = Axis3(
    fig[1,1]; 
    aspect = :data,
    xlabel = L"$r_x$, DU",
    ylabel = L"$r_y$, DU",
    zlabel = L"$r_z$, DU",   
)

# Plot every 10 orbits
for i in 1:10:length(stable_orbits)
    traj = get_full_orbit(stable_orbits[i])
    lines!(ax, traj[1,:], traj[2,:], traj[3,:]; color = :blue)
end

#  Plot the moon
import FileIO
plot_moon(ax, orbit_set) = begin # Defining a function that we can use later
    moon = Sphere(
        Point3f(1.0 - orbit_set.system.mass_ratio, 0.0, 0.0),
        orbit_set.system.radius_secondary / orbit_set.system.DU
    )
    mesh!(ax, moon; color = FileIO.load(CairoMakie.assetpath("moon.png")))
end
plot_moon(ax, stable_orbits)
save("jpl_halos.svg", fig); nothing # hide
```

![](jpl_halos.svg)

# Quick examples
### Plotting a subset of the Earth-Moon Butterfly orbits 
```@example jpl_api
butterfly_orbits = get_jpl_orbits(;
    sys         = "earth-moon",
    family      = "butterfly",
    periodmin   = 12.0,
    periodmax   = 15.0,
    periodunits = "d",
)

butterfly_fig = Figure()
butterfly_ax  = Axis3(
    butterfly_fig[1,1]; 
    aspect = :data,
    xlabel = L"$r_x$, DU",
    ylabel = L"$r_y$, DU",
    zlabel = L"$r_z$, DU",   
)

# Plot every 100 orbits and moon
for i in 1:100:length(butterfly_orbits)
    traj = get_full_orbit(butterfly_orbits[i])
    lines!(butterfly_ax, traj[1,:], traj[2,:], traj[3,:]; color = :blue)
end
plot_moon(butterfly_ax, butterfly_orbits)
save("jpl_butterflys.svg", butterfly_fig) # hide
```

![](jpl_butterflys.svg)

### Plotting all of the Earth-Moon Axial orbits about ``L_4``
```@example jpl_api
axial_orbits = get_jpl_orbits(;
    sys         = "earth-moon",
    family      = "axial",
    libr        = 4,
)

axial_fig = Figure(; size = (800,600))
axial_ax_xy  = Axis(
    axial_fig[1,1]; 
    aspect = DataAspect(),
    xlabel = L"$r_x$, DU",
    ylabel = L"$r_y$, DU",
)
axial_ax_xz  = Axis(
    axial_fig[1,2]; 
    aspect = DataAspect(),
    xlabel = L"$r_x$, DU",
    ylabel = L"$r_z$, DU",
)

# Plot every 200 orbits and moon
for i in 1:200:length(axial_orbits)
    traj = get_full_orbit(axial_orbits[i])
    lines!(axial_ax_xy, traj[1,:], traj[2,:]; color = :blue)
    lines!(axial_ax_xz, traj[2,:], traj[3,:]; color = :blue)
end
save("jpl_axial.svg", axial_fig); nothing # hide
```

![](jpl_axial.svg)

