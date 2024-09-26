"""
    AbstractIndependantVariable

An abstract type to encompas all supported independant variables.
"""
abstract type AbstractIndependantVariable end

"""
    Time

A concrete type to represent time as the independant variable.
"""
struct Time <: AbstractIndependantVariable end

"""
    ArcLength

A concrete type to represent arc-length as the independant variable.
"""
struct ArcLength <: AbstractIndependantVariable end
