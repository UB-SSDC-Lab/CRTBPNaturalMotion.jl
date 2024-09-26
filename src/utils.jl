
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
