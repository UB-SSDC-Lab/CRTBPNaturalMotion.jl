"""
    AbstractInterpolationWrapper

And abstract type that is used to define the interface for interpolation objects.
"""
abstract type AbstractInterpolationWrapper{N} end

"""
    FastChebInterpolation(it::IT) <: AbstractInterpolationWrapper

A wrapper for `FastChebInterp.ChebPoly` that's functional for both the FastChebInterp.jl
implementations of interpolation and least-squares regression (see [FastChebInterp.jl](https://github.com/JuliaMath/FastChebInterp.jl)
for details).

# Fields
- `it::IT`: The `FastChebInterp.ChebPoly` struct.
"""
struct FastChebInterpolation{N, T, Td} <: AbstractInterpolationWrapper{N}
    it::FastChebInterp.ChebPoly{N, T, Td}

    function FastChebInterpolation(
        it::FastChebInterp.ChebPoly{N,T,Td}
    ) where {N, T, Td}
        new{N, T, Td}(it)
    end
end

"""
    save_interp(file::String, it::FastChebInterpolation{IT}) where IT

Save the interpolation object to a file.
"""
function save_interp(file::String, it::FastChebInterpolation)
    save(file, Dict("it" => it.it))
    return nothing
end

"""
    load_interp(file::String) -> FastChebInterpolation{IT}

Load an interpolation object from a file.
"""
function load_interp(file::String)
    it = load(file, "it")
    return FastChebInterpolation(it)
end

"""
    value(m::FastChebInterpolation{N}, τ1[, τ2]) -> SVector{6, Float64}

Evaluate the interpolation at the given `τ1` and `τ2` (if employing a bi-variate
interpolation).

# Arguments
- `m::FastChebInterpolation{IT}`: The interpolation object.
- `τ1`: The variable.
- `τ2`: The second variable (if employing a bi-variate interpolation).

# Returns
- `SVector{6, Float64}`: The interpolated state.
"""
function value(m::FastChebInterpolation{1}, τ1)
    return m.it(τ1)
end
function value(m::FastChebInterpolation{2}, τ1, τ2)
    return m.it(SA[τ1, τ2])
end

"""
    jacobian(m::FastChebInterpolation{1}, τ1) where IT -> SMatrix{6, 1, Float64, 6}

Evaluate the Jacobian of the interpolation at the given `τ1` for a uni-variate interpolation.

# Arguments
- `m::FastChebInterpolation{IT}`: The interpolation object.
- `τ1`: The variable.

# Returns
- `SMatrix{6, 1, Float64, 6}`: The Jacobian of the interpolation.
"""
function jacobian(m::FastChebInterpolation{1}, τ1)
    return chebjacobian(m.it, τ1)[2]
end

"""
    derivative(m::FastChebInterpolation{1}, τ1) where IT -> SVector{6, Float64}

Evaluate the derivative of the interpolation at the given `τ1` for a uni-variate interpolation.

# Arguments
- `m::FastChebInterpolation{IT}`: The interpolation object.
- `τ1`: The variable.

# Returns
- `SVector{6, Float64}`: The derivative of the interpolation.
"""
derivative(m::FastChebInterpolation{1}, τ1) = jacobian(m, τ1)[:]

"""
    jacobian(m::FastChebInterpolation{2}, τ1, τ2) where IT -> SMatrix{6, 2, Float64, 12}

Evaluate the Jacobian of the interpolation at the given `τ1` and `τ2` for a bi-variate interpolation.

# Arguments
- `m::FastChebInterpolation{IT}`: The interpolation object.
- `τ1`: The first variable.
- `τ2`: The second variable.

# Returns
- `SMatrix{6, 2, Float64, 12}`: The Jacobian of the interpolation.
"""
function jacobian(m::FastChebInterpolation{2}, τ1, τ2)
    return chebjacobian(m.it, SA[τ1, τ2])[2]
end
