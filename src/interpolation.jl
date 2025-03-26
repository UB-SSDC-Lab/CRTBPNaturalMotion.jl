"""
    AbstractInterpolationWrapper

And abstract type that is used to define the interface for interpolation objects.
"""
abstract type AbstractInterpolationWrapper{N} end

"""
    INTERP_AD_BACKEND

The backend used for all interpolation AD.
"""
const INTERP_AD_BACKEND = AD.AutoEnzyme(;
    function_annotation=Enzyme.Const, mode=Enzyme.EnzymeCore.Forward
)

"""
    FastChebInterpolation(it::IT) <: AbstractInterpolationWrapper

A wrapper for `FastChebInterp.ChebPoly` that's functional for both the FastChebInterp.jl
implementations of interpolation and least-squares regression (see [FastChebInterp.jl](https://github.com/JuliaMath/FastChebInterp.jl)
for details).

# Fields
- `it::IT`: The `FastChebInterp.ChebPoly` struct.
"""
struct FastChebInterpolation{N,T,Td,P,F<:Function} <: AbstractInterpolationWrapper{N}
    # The chebyshev polynomial
    it::FastChebInterp.ChebPoly{N,T,Td}

    # The differentiated function for computing second derivatives
    fun::F

    # The AD backend for computing second derivatives
    prep::P

    function FastChebInterpolation(it::FastChebInterp.ChebPoly{1,T,Td};) where {T,Td}

        # Define the function and AD preparation
        fun = let it = it
            τ -> derivative(it, τ)
        end
        prep = AD.prepare_derivative(fun, INTERP_AD_BACKEND, zero(Float64))

        return new{1,T,Td,typeof(prep),typeof(fun)}(it, fun, prep)
    end
    function FastChebInterpolation(it::FastChebInterp.ChebPoly{2,T,Td};) where {T,Td}

        # Define the function and AD preparation
        fun = let it = it
            τ -> jacobian(it, τ)
        end
        prep = AD.prepare_jacobian(fun, INTERP_AD_BACKEND, zero(SVector{2,Float64}))

        return new{2,T,Td,typeof(prep),typeof(fun)}(it, fun, prep)
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

function derivative(it::FastChebInterp.ChebPoly{1}, τ1)
    return chebjacobian(it, τ1)[2][:]
end
function jacobian(it::FastChebInterp.ChebPoly, τ)
    return chebjacobian(it, τ)[2]
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
    return jacobian(m.it, τ1)
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
derivative(m::FastChebInterpolation{1}, τ1) = derivative(m.it, τ1)

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
    return jacobian(m.it, SA[τ1, τ2])
end

"""
    second_derivative(m::FastChebInterpolation{1}, τ1) -> SVector{6, Float64}

Evaluate the second derivative of the interpolation at the given `τ1` for a
uni-variate interpolation. The second derivatives are computed by employing Enzyme to
differentiate the analytical computation of the first derivative.

# Arguments
- `m::FastChebInterpolation{IT}`: The interpolation object.
- `τ1`: The variable.

# Returns
- `SVector{6, Float64}`: The second derivative of the interpolation.
"""
function second_derivative(m::FastChebInterpolation{1}, τ1)
    return AD.derivative(m.fun, m.prep, INTERP_AD_BACKEND, τ1)
end

"""
    hessian(m::FastChebInterpolation{1}, τ1, τ2) -> SMatrix{6, 1, Float64, 6}

Evaluate the Hessian of the interpolation at the given `τ1`. This computes the second
derivative of the interpolation at the given `τ1` for a uni-variate interpolation and returns
a SMatrix{6, 1, Float64, 6} (rather than an SVector{6, Float64} like `second_derivative`).

# Arguments
- `m::FastChebInterpolation{IT}`: The interpolation object.
- `τ1`: The variable.

# Returns
- `SMatrix{6, 1, Float64, 6}`: The Hessian of the interpolation.
"""
hessian(m::FastChebInterpolation{1}, τ1) = SMatrix{6,1,Float64,6}(second_derivative(m, τ1))

"""
    hessian(m::FastChebInterpolation{2}, τ1, τ2) -> SArray{Tuple{6, 2, 2}, Float64, 3, 24}

Evaluate the Hessian of the interpolation at the given `τ1` and `τ2`. This computes
the second derivative of the interpolation at the given `τ1` and `τ2` for a bi-variate
interpolation and returns a SArray{Tuple{6, 2, 2}, Float64, 3, 24}.

# Arguments
- `m::FastChebInterpolation{IT}`: The interpolation object.
- `τ1`: The first variable.
- `τ2`: The second variable.

# Returns
- `SArray{Tuple{6, 2, 2}, Float64, 3, 24}`: The Hessian of the interpolation.
"""
function hessian(m::FastChebInterpolation{2}, τ1, τ2)
    hess_reshaped = AD.jacobian(m.fun, m.prep, INTERP_AD_BACKEND, SA[τ1, τ2])
    return hess_reshaped.parent
end
