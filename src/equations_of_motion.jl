
"""
    natural_crtbp_eom(x, p, t)

Compute the equations of motion for the CRTBP model with only natural motion.

# Arguments
- `x::SVector{6}`: State vector, i.e., [position^T; velocity^T].
- `p::Tuple{Float64,Float64,Float64}`: Tuple of parameters, (mu,)
- `t::Float64`: Time.
"""
function natural_crtbp_eom(x, p, t)
    # Grab states 
    rx = x[1]
    ry = x[2]
    rz = x[3]
    vx = x[4]
    vy = x[5]
    vz = x[6]

    # Grab parameters
    mu = p[1]

    # Start of matlab generated code
    # See ./code_gen/crtbp_costate_eom_code_gen.m for details
    t2 = mu + rx
    t3 = ry * ry
    t4 = rz * rz
    t6 = mu - 1.0
    t7 = t2 - 1.0
    t8 = t2 * t2
    t9 = t7 * t7
    t10 = t3 + t4 + t8
    t11 = t3 + t4 + t9

    #t12 = 1.0/pow(t10,3.0/2.0);
    tt1 = sqrt(t10)
    tt13 = tt1 * tt1 * tt1
    t12 = 1.0 / tt13

    #t13 = 1.0/pow(t11,3.0/2.0);
    tt2 = sqrt(t11)
    tt23 = tt2 * tt2 * tt2
    t13 = 1.0 / tt23

    # Compute and return dynamics
    return SA[
        vx,
        vy,
        vz,
        rx + vy * 2.0 - mu * t7 * t13 + t2 * t6 * t12,
        ry - vx * 2.0 - mu * ry * t13 + ry * t6 * t12,
        -mu * rz * t13 + rz * t6 * t12,
    ]
end

"""
    natural_crtbp_eom_with_stm(x, p, t)

Compute the equations of motion for the CRTBP model with the state transition matrix.

# Arguments
- `x::SMatrix{6,7,Float64}`: Matrix with first column ast the state and the remaining 6 x 6 as the STM.
- `p::Tuple{Float64,}`: Tuple of parameters, (mu,)
- `t::Float64`: Time.
"""
function natural_crtbp_eom_with_stm(z, p, t)
    # Grab states
    rx = z[1, 1]
    ry = z[2, 1]
    rz = z[3, 1]
    vx = z[4, 1]
    vy = z[5, 1]
    vz = z[6, 1]

    # Grab STM
    STM = z[SA[1, 2, 3, 4, 5, 6], SA[2, 3, 4, 5, 6, 7]]

    # Grab parameters
    mu = p[1]

    # Start of matlab generated code
    t2 = mu + rx
    t5 = ry * ry
    t6 = rz * rz
    t7 = mu - 1.0
    t8 = t2 - 1.0
    t9 = t2 * t2
    t10 = t8 * t8
    t13 = t5 + t6 + t9
    t14 = t5 + t6 + t10

    #t15 = 1.0/pow(t13,3.0/2.0);
    #t16 = 1.0/pow(t13,5.0/2.0);
    tt1 = sqrt(t13)
    tt13 = tt1 * tt1 * tt1
    tt15 = tt13 * tt1 * tt1
    t15 = 1.0 / tt13
    t16 = 1.0 / tt15

    #t17 = 1.0/pow(t14,3.0/2.0);
    #t18 = 1.0/pow(t14,5.0/2.0);
    tt2 = sqrt(t14)
    tt23 = tt2 * tt2 * tt2
    tt25 = tt23 * tt2 * tt2
    t17 = 1.0 / tt23
    t18 = 1.0 / tt25

    t20 = t7 * t15
    t23 = ry * rz * t7 * t16 * 3.0
    t19 = mu * t17
    t22 = mu * ry * rz * t18 * 3.0
    t24 = -t23
    t21 = -t19
    t25 = t22 + t24

    # State dynamics
    dx = SA[
        vx,
        vy,
        vz,
        rx + vy * 2.0 + t2 * t20 + t8 * t21,
        ry - vx * 2.0 + ry * t20 + ry * t21,
        rz * t20 + rz * t21,
    ]

    # Dynamics Jacobian
    F14 = 1.0
    F25 = 1.0
    F36 = 1.0
    F41 = t20 + t21 + mu * t10 * t18 * 3.0 - t7 * t9 * t16 * 3.0 + 1.0
    F42 = mu * ry * t8 * t18 * 3.0 - ry * t2 * t7 * t16 * 3.0
    F43 = mu * rz * t8 * t18 * 3.0 - rz * t2 * t7 * t16 * 3.0
    F45 = 2.0
    F51 = mu * ry * t8 * t18 * 3.0 - ry * t2 * t7 * t16 * 3.0
    F52 = t20 + t21 + mu * t5 * t18 * 3.0 - t5 * t7 * t16 * 3.0 + 1.0
    F53 = t25
    F54 = -2.0
    F61 = mu * rz * t8 * t18 * 3.0 - rz * t2 * t7 * t16 * 3.0
    F62 = t25
    F63 = t20 + t21 + mu * t6 * t18 * 3.0 - t6 * t7 * t16 * 3.0
    F = SA[
        0.0 0.0 0.0 F14 0.0 0.0
        0.0 0.0 0.0 0.0 F25 0.0
        0.0 0.0 0.0 0.0 0.0 F36
        F41 F42 F43 0.0 F45 0.0
        F51 F52 F53 F54 0.0 0.0
        F61 F62 F63 0.0 0.0 0.0
    ]

    # Compute stm dynamics
    dSTM = F * STM

    # Return 
    return SMatrix{6,7}(
        dx[1],
        dx[2],
        dx[3],
        dx[4],
        dx[5],
        dx[6],
        dSTM[1],
        dSTM[2],
        dSTM[3],
        dSTM[4],
        dSTM[5],
        dSTM[6],
        dSTM[7],
        dSTM[8],
        dSTM[9],
        dSTM[10],
        dSTM[11],
        dSTM[12],
        dSTM[13],
        dSTM[14],
        dSTM[15],
        dSTM[16],
        dSTM[17],
        dSTM[18],
        dSTM[19],
        dSTM[20],
        dSTM[21],
        dSTM[22],
        dSTM[23],
        dSTM[24],
        dSTM[25],
        dSTM[26],
        dSTM[27],
        dSTM[28],
        dSTM[29],
        dSTM[30],
        dSTM[31],
        dSTM[32],
        dSTM[33],
        dSTM[34],
        dSTM[35],
        dSTM[36],
    )
end

"""
    natural_crtbp_eom_with_arclen(x, p, t)

Compute the equations of motion for the CRTBP model with only natural motion, 
employing the arclength along the trajectory as an additional state.

# Arguments
- `x::SVector{7}`: State vector, i.e., [position^T; velocity^T; arc-len].
- `p::Tuple{Float64,}`: Tuple of parameters, (mu,)
- `t::Float64`: Time.
"""
function natural_crtbp_eom_with_arclen(x, p, t)
    # Grab states 
    rx = x[1]
    ry = x[2]
    rz = x[3]
    vx = x[4]
    vy = x[5]
    vz = x[6]

    # Grab parameters
    mu = p[1]

    # Start of matlab generated code
    # See ./code_gen/crtbp_costate_eom_code_gen.m for details
    t2 = mu + rx
    t3 = ry * ry
    t4 = rz * rz
    t6 = mu - 1.0
    t7 = t2 - 1.0
    t8 = t2 * t2
    t9 = t7 * t7
    t10 = t3 + t4 + t8
    t11 = t3 + t4 + t9

    #t12 = 1.0/pow(t10,3.0/2.0);
    tt1 = sqrt(t10)
    tt13 = tt1 * tt1 * tt1
    t12 = 1.0 / tt13

    #t13 = 1.0/pow(t11,3.0/2.0);
    tt2 = sqrt(t11)
    tt23 = tt2 * tt2 * tt2
    t13 = 1.0 / tt23

    # Compute and return dynamics
    return SA[
        vx,
        vy,
        vz,
        rx + vy * 2.0 - mu * t7 * t13 + t2 * t6 * t12,
        ry - vx * 2.0 - mu * ry * t13 + ry * t6 * t12,
        -mu * rz * t13 + rz * t6 * t12,
        sqrt(vx * vx + vy * vy + vz * vz),
    ]
end

"""
    natural_crtbp_eom_with_independant_arclen(x, p, t)

Compute the equations of motion for the CRTBP model with only natural motion
where the independant variable has been transformed from time to arc-length.

# Arguments
- `x::SVector{7}`: State vector, i.e., [position^T; velocity^T].
- `p::Tuple{Float64,}`: Tuple of parameters, (mu,)
- `t::Float64`: Time.
"""
function natural_crtbp_eom_with_independant_arclen(x, p, t)
    # Grab states 
    rx = x[1]
    ry = x[2]
    rz = x[3]
    vx = x[4]
    vy = x[5]
    vz = x[6]

    # Grab parameters
    mu = p[1]

    # Start of matlab generated code
    # See ./code_gen/crtbp_costate_eom_code_gen.m for details
    t2 = mu + rx
    t3 = ry * ry
    t4 = rz * rz
    t5 = vx * vx
    t6 = vy * vy
    t7 = vz * vz
    t8 = mu - 1.0
    t9 = t2 - 1.0
    t10 = t2 * t2
    t12 = t5 + t6 + t7
    t11 = t9 * t9
    t13 = t3 + t4 + t10
    t15 = 1.0 / sqrt(t12)
    t14 = t3 + t4 + t11

    #t16 = 1.0/pow(t13,3.0/2.0);
    tt1 = sqrt(t13)
    tt13 = tt1 * tt1 * tt1
    t16 = 1.0 / tt13

    #t17 = 1.0/pow(t14,3.0/2.0);
    tt2 = sqrt(t14)
    tt23 = tt2 * tt2 * tt2
    t17 = 1.0 / tt23

    return SA[
        t15 * vx,
        t15 * vy,
        t15 * vz,
        t15 * (rx + vy * 2.0 - mu * t9 * t17 + t2 * t8 * t16),
        t15 * (ry - vx * 2.0 - mu * ry * t17 + ry * t8 * t16),
        -t15 * (mu * rz * t17 - rz * t8 * t16),
    ]
end

"""
    natural_crtbp_and_time_eom_with_independant_arclen(x, p, t)

Compute the equations of motion for the CRTBP model with an extra state 
for the time with only natural motion. Here the independant variable has been 
transformed from time to arc-length.

# Arguments
- `x::SVector{7}`: State vector, i.e., [position^T; velocity^T].
- `p::Tuple{Float64,}`: Tuple of parameters, (mu,)
- `t::Float64`: Time.
"""
function natural_crtbp_and_time_eom_with_independant_arclen(x, p, t)
    # Grab states 
    rx = x[1]
    ry = x[2]
    rz = x[3]
    vx = x[4]
    vy = x[5]
    vz = x[6]

    # Grab parameters
    mu = p[1]

    # Start of matlab generated code
    # See ./code_gen/crtbp_costate_eom_code_gen.m for details
    t2 = mu + rx
    t3 = ry * ry
    t4 = rz * rz
    t5 = vx * vx
    t6 = vy * vy
    t7 = vz * vz
    t8 = mu - 1.0
    t9 = t2 - 1.0
    t10 = t2 * t2
    t12 = t5 + t6 + t7
    t11 = t9 * t9
    t13 = t3 + t4 + t10
    t15 = 1.0 / sqrt(t12)
    t14 = t3 + t4 + t11

    #t16 = 1.0/pow(t13,3.0/2.0);
    tt1 = sqrt(t13)
    tt13 = tt1 * tt1 * tt1
    t16 = 1.0 / tt13

    #t17 = 1.0/pow(t14,3.0/2.0);
    tt2 = sqrt(t14)
    tt23 = tt2 * tt2 * tt2
    t17 = 1.0 / tt23

    return SA[
        t15 * vx,
        t15 * vy,
        t15 * vz,
        t15 * (rx + vy * 2.0 - mu * t9 * t17 + t2 * t8 * t16),
        t15 * (ry - vx * 2.0 - mu * ry * t17 + ry * t8 * t16),
        -t15 * (mu * rz * t17 - rz * t8 * t16),
        t15,
    ]
end
