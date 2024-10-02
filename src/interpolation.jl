"""
    AbstractInterpolationWrapper

And abstract type that is used to define the interface for interpolation objects.
"""
abstract type AbstractInterpolationWrapper end

"""
    FastChebInterpolation(it::IT) <: AbstractInterpolationWrapper

A wrapper for `FastChebInterp.ChebPoly` that's functional for both the FastChebInterp.jl
implementations of interpolation and least-squares regression (see [FastChebInterp.jl](https://github.com/JuliaMath/FastChebInterp.jl)
for details).

# Fields
- `it::IT`: The `FastChebInterp.ChebPoly` struct.
"""
struct FastChebInterpolation{IT <: FastChebInterp.ChebPoly} <: AbstractInterpolationWrapper
    it::IT

    function FastChebInterpolation(it::IT) where IT
        new{IT}(it)
    end
end

"""
    value(m::FastChebInterpolation{IT}, τ1[, τ2]) where IT -> SVector{6, Float64}

Evaluate the interpolation at the given `τ1` and `τ2` (if employing a bi-variate
interpolation).

# Arguments
- `m::FastChebInterpolation{IT}`: The interpolation object.
- `τ1`: The variable.
- `τ2`: The second variable (if employing a bi-variate interpolation).

# Returns
- `SVector{6, Float64}`: The interpolated state.
"""
function value(m::FastChebInterpolation{IT}, τ1) where IT
    return m.it(τ1)
end
function value(m::FastChebInterpolation{IT}, τ1, τ2) where IT
    return m.it(SA[τ1, τ2])
end

"""
    jacobian(m::FastChebInterpolation{IT}, τ1) where IT -> SVector{6, Float64}

Evaluate the Jacobian of the interpolation at the given `τ1` for a uni-variate interpolation.

# Arguments
- `m::FastChebInterpolation{IT}`: The interpolation object.
- `τ1`: The variable.

# Returns
- `SVector{6, Float64}`: The Jacobian of the interpolation.
"""
function jacobian(m::FastChebInterpolation{IT}, τ1) where IT
    return chebjaconian(m.it, τ1)[2]
end

"""
    jacobian(m::FastChebInterpolation{IT}, τ1, τ2) where IT -> SMatrix{6, 2, Float64, 12}

Evaluate the Jacobian of the interpolation at the given `τ1` and `τ2` for a bi-variate interpolation.

# Arguments
- `m::FastChebInterpolation{IT}`: The interpolation object.
- `τ1`: The first variable.
- `τ2`: The second variable.

# Returns
- `SMatrix{6, 2, Float64, 12}`: The Jacobian of the interpolation.
"""
function jacobian(m::FastChebInterpolation{IT}, τ1, τ2) where IT
    return chebjaconian(m.it, SA[τ1, τ2])[2]
end
