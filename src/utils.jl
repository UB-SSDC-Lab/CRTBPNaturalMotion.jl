
"""
    initial_state_with_stm(x0::SVector{6,T}) where T

Create the initial state with the state transition matrix (STM).

# Arguments
- `x0::SVector{6,T}`: Initial state vector.

# Returns
- `SMatrix{6,7,T,42}`: Matrix with first column as the state and the remaining 6 x 6 as
    the STM initial condition (i.e., [x0, I]).
"""
@inline function initial_state_with_stm(x0::SVector{6,T}) where T
    return SMatrix{6,7,T,42}(
        x0[1], x0[2], x0[3], x0[4], x0[5], x0[6],
        1.0,  0.0, 0.0,  0.0, 0.0,  0.0,
        0.0,  1.0, 0.0,  0.0, 0.0,  0.0,
        0.0,  0.0, 1.0,  0.0, 0.0,  0.0,
        0.0,  0.0, 0.0,  1.0, 0.0,  0.0,
        0.0,  0.0, 0.0,  0.0, 1.0,  0.0,
        0.0,  0.0, 0.0,  0.0, 0.0,  1.0,
    )
end

"""
    get_stm(z::SMatrix{6,7,T,42}) where T

Get the state transition matrix (STM) from the full integration matrix.

# Arguments
- `z::SMatrix{6,7,T,42}`: Matrix with first column as the state and the remaining 6 x 6 as
    the STM.

# Returns
- `SMatrix{6,6,T,36}`: The state transition matrix.
"""
@inline function get_stm(z::SMatrix{6,7,T,42}) where {T}
    return z[SA[1,2,3,4,5,6],SA[2,3,4,5,6,7]]
end

"""
    get_state_and_stm(z::SMatrix{6,7,T,42}) where T

Get the state and state transition matrix (STM) from the full integration matrix.

# Arguments
- `z::SMatrix{6,7,T,42}`: Matrix with first column as the state and the remaining 6 x 6 as
    the STM.

# Returns
- `Tuple{SVector{6,T},SMatrix{6,6,T,36}`: The state vector and the state transition matrix
    in a return tuple.
"""
@inline function get_state_and_stm(z::SMatrix{6,7,T,42}) where {T}
    return z[SA[1,2,3,4,5,6]], get_stm(z)
end

"""
    jacobi_integral(x, μ)

Computes the Jacobi integral for the CRTBP given the state `x` and mass parameter `μ`.

# Arguments
- `x::AbstractVector`: State vector.
- `μ::Real`: Mass parameter.
"""
function jacobi_integral(x, μ)
    xr1 = x[1] + μ
    xr2 = xr1 - 1.0
    r1 = sqrt(xr1*xr1 + x[2]*x[2] + x[3]*x[3])
    r2 = sqrt(xr2*xr2 + x[2]*x[2] + x[3]*x[3])
    return (x[1]*x[1] + x[2]*x[2]) + 2.0*((1.0 - μ)/r1 + μ/r2) - (x[4]*x[4] + x[5]*x[5] + x[6]*x[6])
end

"""
    jacobi_integral_gradient(x, μ)

Computes the gradient of the Jacobi integral for the CRTBP given the state `x` and mass
parameter `μ`.

# Arguments
- `x::AbstractVector`: State vector.
- `μ::Real`: Mass parameter.
"""
function jacobi_integral_gradient(x, mu)
    # Grab states
    rx = x[1]; ry = x[2]; rz = x[3]
    vx = x[4]; vy = x[5]; vz = x[6]

    # Matlab generated code
    t2 = mu+rx;
    t4 = rx*2.0;
    t5 = ry*ry;
    t6 = rz*rz;
    t7 = mu-1.0;
    t8 = t2-1.0;
    t9 = t2*t2;
    t10 = t8*t8;
    t11 = t5+t6+t9;
    t12 = t5+t6+t10;

    #t13 = 1.0/pow(t11,3.0/2.0);
    tt1 = sqrt(t11)
    tt13 = tt1*tt1*tt1
    t13 = 1.0/tt13

    #t14 = 1.0/pow(t12,3.0/2.0);
    tt2 = sqrt(t12)
    tt23 = tt2*tt2*tt2
    t14 = 1.0/tt23

    return SA[
        t4-mu*t8*t14*2.0+t2*t7*t13*2.0,
        ry*2.0-mu*ry*t14*2.0+ry*t7*t13*2.0,
        mu*rz*t14*-2.0+rz*t7*t13*2.0,
        vx*-2.0,
        vy*-2.0,
        vz*-2.0,
    ]
end
