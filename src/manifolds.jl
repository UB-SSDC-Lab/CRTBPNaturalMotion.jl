abstract type AbstractManifold end

struct InvariantManifold{PO <: AbstractPeriodicOrbit} <: AbstractManifold
    # The periodic orbit
    orbit::PO

    # The manifold perturbation size (change in position)
end
