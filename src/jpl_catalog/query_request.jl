
function request_data(query::String)
    # HTTP request
    r = HTTP.request("GET", query)

    # Check status code
    if r.status != 200
        if r.status == 400
            error("Bad Request: The request constained invalid keywords and/or content.")
        elseif r.status == 405
            error("Method Not Allowed: The request used an incorrect method.")
        elseif r.status == 500
            error(
                "Internal Server Error: The database is not available at the time of the request.",
            )
        elseif r.status == 503
            error(
                "Service Unavailable: The server is currently unable to handle the request due to a temporary overloading or maintenance of the server.",
            )
        else
            error("Unknown HTTP status code: $(r.status)")
        end
    end

    # Parse JSON string and return
    return JSON3.read(String(r.body), JSONOrbitData; parsequoted=true)
end
