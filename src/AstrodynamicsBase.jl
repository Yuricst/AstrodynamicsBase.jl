module AstrodynamicsBase

    using LinearAlgebra
    using SPICE
    using Dates

    include("misc.jl")
    include("trigonometry.jl")
    include("transformations.jl")
    include("body_properties.jl")
    include("elements.jl")
    include("canonical.jl")
    include("spice_helper.jl")
    include("canonical_thrust.jl")
    include("kepler_problem.jl")
    include("diffeq_twobody.jl")
    include("orbit.jl")
    include("tle.jl")

    export mod_custom, angle_difference
    export acos_safe, asin_safe
    export _rotmat_ax1,
        _rotmat_ax2,
        _rotmat_ax3,
        rotating2inertial, inertial2rotating, perifocal2geocentric
    export get_gm,
        get_gm_de431,
        get_body_radius,
        get_body_period,
        get_body_sma,
        get_body_soi,
        get_canonical_param,
        get_semiMajorAxes
    export kep2cart, cart2kep, kep2mee, mee2kep, mee2cart, cart2mee, cart2poincare, get_period
    export spice_furnsh,
        spice_utc2et,
        spice_et2utc,
        et2datetime,
        get_spice_locate_function,
        interp_spkssb,
        shift_eclipj2000_to_bodycenter,
        spice_spkssb
    export dimensional2canonical_thrust,
        dimensional2canonical_mdot,
        dimensional2canonical_isp
    export tle_from_file

end # module
