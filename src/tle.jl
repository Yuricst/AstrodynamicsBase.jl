"""
Handling TLEs
"""

"""
    tle_three_lines(sat_string::Vector{String})

Process vector of three strings corresponding to Name, Line 1, Line 2 of TLE file.
"""
function tle_three_lines(sat_string::Vector{String})
    # proxy for first and second lines
    line1 = sat_string[2]
    line2 = sat_string[3]

    dict_per_sat = Dict(
        "Name" => strip(sat_string[1]),
        # first line
        "CatalogNumber"           => parse(Int64, line1[3:7]),
        "Clasification"           => line1[8],
        "InternationalDesignator" => strip(line1[10:17]),
        "EpochYear"               => parse(Int64, line1[19:20]),
        "EopchDay"                => parse(Float64, line1[21:32]),
        "dMeanMotion"             => parse(Float64, line1[34:43]),
        "d2MeanMotion"            => line1[45:52],      # FIXME
        "BSTAR"                   => line1[54:61],      # FIXME
        "EphemerisType"           => parse(Int64, line1[63]),
        "ElementSetNumber"        => parse(Int64, line1[65:68]),
        # second line
        "Inclination"       => parse(Float64, line2[9:16]),
        "RAAN"              => parse(Float64, line2[18:25]),
        "Eccentricity"      => parse(Float64, "0."*line2[27:33]),
        "ArgumentOfPerigee" => parse(Float64, line2[35:42]),
        "MeanAnomaly"       => parse(Float64, line2[44:51]),
        "MeanMotion"        => parse(Float64, line2[53:63]),
        "RevNumberAtEpoch"  => parse(Int64, line2[64:68]),
    )
    return dict_per_sat
end


"""
    tle_from_file(filename::String)

Get dictionary of satellites from TLE file, with keys corresponding to their names.
Input text file expected in the format of what can be obtained from celestrak.
"""
function tle_from_file(filename::String)
    # load list of tles
    string_list = readlines(filename)
    # separate into each satellite (3 lines per satellite)
    string_sats = [string_list[1+3*(i-1):3*i] for i = 1:div(length(string_list),3)];
    # for each satellite, get dictionary
    sat_dict_list = [
        tle_three_lines(string_sats[k]) for k = 1:length(string_sats)
    ]
    # create dictionary by name
    sat_dict = Dict()
    for sat in sat_dict_list
        sat_dict[sat["Name"]] = sat
    end
    return sat_dict
end


"""
    tle2kep(tle_dict::Dict, μ_Earth::Real=398600.4418, day::Real=86400.002)

Get Keplerian elements from TLE dictionary entry.
Note the definition of a day in TLE is mean solar day, not sidereal day.
"""
function tle2kep(tle_dict::Dict, μ_Earth::Real=398600.4418, day::Real=86400.002)
    n = tle_dict["MeanMotion"]*2π/day  # rev/day -> rad/sec
    sma = (μ_Earth/n^2)^(1/3)
    TA = anomaly_mean_to_true(deg2rad(tle_dict["MeanAnomaly"]), tle_dict["Eccentricity"])
    kep_elts = [
        sma,
        tle_dict["Eccentricity"],
        deg2rad(tle_dict["Inclination"]),
        deg2rad(tle_dict["RAAN"]),
        deg2rad(tle_dict["ArgumentOfPerigee"]),
        TA
    ]
    return kep_elts
end


"""
    query(search_key::String, dictionary::Dict)

Query dictionary based on partial key
"""
function query(search_key::String, dictionary::Dict)
    query_result = Dict()
    for key in keys(dictionary)
        if occursin(search_key, key)
            query_result[key] = dictionary[key]
        end
    end
    return query_result
end
