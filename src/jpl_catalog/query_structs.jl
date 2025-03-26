# ===== Structs for deserializing JSON data using JSON3.jl

# NOTE: These structs have a 1-to-1 correspondence with the JSON data returned
# by the JPL API and should generally not be used by the user

mutable struct JSONSignature
    source::String
    version::String
    JSONSignature() = new("N/A", "N/A")
end

mutable struct JSONSystem
    name::String
    mass_ratio::Float64
    radius_secondary::Float64
    L1::Vector{Float64}
    L2::Vector{Float64}
    L3::Vector{Float64}
    L4::Vector{Float64}
    L5::Vector{Float64}
    lunit::Float64
    tunit::Float64
    function JSONSystem()
        return new(
            "N/A",
            NaN,
            NaN,
            Vector{Float64}(undef, 0),
            Vector{Float64}(undef, 0),
            Vector{Float64}(undef, 0),
            Vector{Float64}(undef, 0),
            Vector{Float64}(undef, 0),
            NaN,
            NaN,
        )
    end
end

mutable struct JSONLimits
    jacobi::Vector{Float64}
    period::Vector{Float64}
    stability::Vector{Float64}
    function JSONLimits()
        return new(
            Vector{Float64}(undef, 0), Vector{Float64}(undef, 0), Vector{Float64}(undef, 0)
        )
    end
end

mutable struct JSONOrbitData
    signature::JSONSignature
    system::JSONSystem
    family::String
    libration_point::Int
    branch::String
    limits::JSONLimits
    count::Int
    fields::Vector{String}
    data::Vector{Vector{Float64}}
    function JSONOrbitData()
        return new(
            JSONSignature(),
            JSONSystem(),
            "N/A",
            0,
            "N/A",
            JSONLimits(),
            0,
            ["x", "y", "z", "vx", "vy", "vz", "jacobi", "period", "stability"],
            Vector{Vector{Float64}}(undef, 0),
        )
    end
end

StructTypes.StructType(::Type{JSONSignature}) = StructTypes.Mutable()
StructTypes.StructType(::Type{JSONSystem}) = StructTypes.Mutable()
StructTypes.StructType(::Type{JSONLimits}) = StructTypes.Mutable()
StructTypes.StructType(::Type{JSONOrbitData}) = StructTypes.Mutable()
