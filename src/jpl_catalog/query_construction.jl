# Define api link
const global api = "https://ssd-api.jpl.nasa.gov/periodic_orbits.api"

# Define allowable query parameters
const global valid_sys = [
    "earth-moon",
    "mars-phobos",
    "sun-earth",
    "sun-mars",
    "jupiter-europa",
    "saturn-enceladus",
    "saturn-titan",
]
const global valid_family = [
    "halo",
    "vertical",
    "axial",
    "lyapunov",
    "longp",
    "short",
    "butterfly",
    "dragonfly",
    "resonant",
    "dro",
    "dpo",
    "lpo",
]
const global valid_periodunits = ["s", "h", "d", "TU"]

function construct_query(;
    sys::String="earth-moon",
    family::String="halo",
    libr::Int=1,
    branch::String="N",
    periodmin::Union{Float64,Nothing}=nothing,
    periodmax::Union{Float64,Nothing}=nothing,
    periodunits::Union{String,Nothing}=nothing,
    jacobimin::Union{Float64,Nothing}=nothing,
    jacobimax::Union{Float64,Nothing}=nothing,
    stabmin::Union{Float64,Nothing}=nothing,
    stabmax::Union{Float64,Nothing}=nothing,
)
    # Check sys and family arguments
    if !(sys in valid_sys)
        # Construct error message
        err_msg = "Invalid sys: $sys.\n"
        err_msg *= "  Options include:\n"
        for (i, sys_option) in enumerate(valid_sys)
            err_msg *= "    - \"$sys_option\""
            if i < length(valid_sys)
                err_msg *= "\n"
            end
        end
        throw(ArgumentError(err_msg))
    end
    if !(family in valid_family)
        err_msg = "Invalid family: $family.\n"
        err_msg *= "  Options include:\n"
        for (i, fam) in enumerate(valid_family)
            err_msg *= "    - \"$fam\""
            if i < length(valid_sys)
                err_msg *= "\n"
            end
        end
        throw(ArgumentError(err_msg))
    end

    # Check periodunits argument
    if periodunits !== nothing
        if !(periodunits in valid_periodunits)
            err_msg = "Invalid periodunits: $periodunits.\n"
            err_msg *= "  Options include:\n"
            for (i, unit) in enumerate(valid_periodunits)
                err_msg *= "    - \"$unit\""
                if i < length(valid_periodunits)
                    err_msg *= "\n"
                end
            end
            throw(ArgumentError(err_msg))
        end
    end

    # Check libration point argument
    if family in ["lyapunov", "halo"]
        if libr < 1 || libr > 3
            throw(ArgumentError("Invalid libration point: $libr. Must be 1, 2, or 3."))
        end
    elseif family in ["longp", "short"]
        if libr < 4 || libr > 5
            throw(ArgumentError("Invalid libration point: $libr. Must be 4 or 5."))
        end
    elseif family in ["axial", "vertical"]
        if libr < 1 || libr > 5
            throw(
                ArgumentError("Invalid libration point: $libr. Must be 1, 2, 3, 4, or 5.")
            )
        end
    end

    # Check branch argument
    if family in ["halo", "dragonfly", "butterfly"]
        if !(branch in ["N", "S"])
            throw(ArgumentError("Invalid branch: $branch. Must be \"N\" or \"S\"."))
        end
    elseif family == "lpo"
        if !(branch in ["E", "W"])
            throw(ArgumentError("Invalid branch: $branch. Must be \"E\" or \"W\"."))
        end
    elseif family == "resonant"
        # This check can probably be improved
        try
            parse(Int, branch)
        catch
            throw(
                ArgumentError(
                    "Invalid branch: $branch. Must be parsable as an integer, i.e., for \"12\" for 1:2.",
                ),
            )
        end
    end

    # Construct and return query
    if family == "halo" # Need libr, branch
        query = api * "?sys=$sys&family=halo&libr=$libr&branch=$branch"
    elseif family == "vertical" # Need libr
        query = api * "?sys=$sys&family=vertical&libr=$libr"
    elseif family == "axial" # Need libr
        query = api * "?sys=$sys&family=axial&libr=$libr"
    elseif family == "lyapunov" # Need libr
        query = api * "?sys=$sys&family=lyapunov&libr=$libr"
    elseif family == "longp" # Need libr
        query = api * "?sys=$sys&family=longp&libr=$libr"
    elseif family == "short" # Need libr
        query = api * "?sys=$sys&family=short&libr=$libr"
    elseif family == "butterfly" # Need branch
        query = api * "?sys=$sys&family=butterfly&branch=$branch"
    elseif family == "dragonfly" # Need branch
        query = api * "?sys=$sys&family=dragonfly&branch=$branch"
    elseif family == "resonant" # Need branch
        query = api * "?sys=$sys&family=resonant&branch=$branch"
    elseif family == "dro"
        query = api * "?sys=$sys&family=dro"
    elseif family == "dpo"
        query = api * "?sys=$sys&family=dpo"
    elseif family == "lpo" # Need branch
        query = api * "?sys=$sys&family=lpo&branch=$branch"
    end

    # Add filters
    if periodmin !== nothing
        query *= "&periodmin=$periodmin"
    end
    if periodmax !== nothing
        query *= "&periodmax=$periodmax"
    end
    if periodunits !== nothing && (periodmin !== nothing || periodmax !== nothing)
        query *= "&periodunits=$periodunits"
    end
    if jacobimin !== nothing
        query *= "&jacobimin=$jacobimin"
    end
    if jacobimax !== nothing
        query *= "&jacobimax=$jacobimax"
    end
    if stabmin !== nothing
        query *= "&stabmin=$stabmin"
    end
    if stabmax !== nothing
        query *= "&stabmax=$stabmax"
    end

    return query
end
