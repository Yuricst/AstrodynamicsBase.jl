module AstrodynamicsBase

    using LinearAlgebra
    using SPICE

    include("trigonometry.jl")
    include("body_properties.jl")
    include("elements.jl")
    include("canonical.jl")
    include("spice_helper.jl")
    include("canonical_thrust.jl")

    export acos_safe, asin_safe
    export get_body_soi, get_canonical_param
    export kep2cart, cart2kep, kep2mee, mee2kep, mee2cart, cart2mee, cart2poincare
    export spice_furnsh,
        spice_utc2et,
        spice_et2utc,
        get_spice_locate_function,
        interp_spkssb,
        shift_eclipj2000_to_bodycenter,
        spice_spkssb
    export dimensional2canonical_thrust, 
        dimensional2canonical_mdot, 
        dimensional2canonical_isp

end # module
