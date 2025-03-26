# Function to get the Earth-Moon mass parameter
get_earth_moon_mass_parameter() = 1.21506038e-2

# Functions to get scenario orbit initial states
get_L1_halo_ic() = SA[0.8233851820, 0.0, -0.0222775563, 0.0, 0.1341841703, 0.0]
function get_L1aziz_halo_ic()
    return SA[0.84842330082624, 0.0, 0.17351888331464177, 0.0, 0.2636116677034408, 0.0]
end
function get_L2aziz_halo_ic()
    return SA[1.1601775563986034, 0.0, -0.12427171585626111, 0.0, -0.20845667597559503, 0.0]
end

function generate_L1_halo()
    cond(x, t, integ) = x[2]
    affect!(integ) = integ.u[3] < 0.0 && terminate!(integ)
    return propagate_return_all_states(
        get_L1_halo_ic(), (0.0, Inf), get_earth_moon_mass_parameter(), cond, affect!, Time;
    )
end

function generate_L1aziz_halo()
    cond(x, t, integ) = x[2]
    affect!(integ) = integ.u[3] > 0.0 && terminate!(integ)
    return propagate_return_all_states(
        get_L1aziz_halo_ic(),
        (0.0, Inf),
        get_earth_moon_mass_parameter(),
        cond,
        affect!,
        Time;
    )
end

function generate_L2aziz_halo()
    cond(x, t, integ) = x[2]
    affect!(integ) = integ.u[3] < 0.0 && terminate!(integ)
    return propagate_return_all_states(
        get_L2aziz_halo_ic(),
        (0.0, Inf),
        get_earth_moon_mass_parameter(),
        cond,
        affect!,
        Time;
    )
end
