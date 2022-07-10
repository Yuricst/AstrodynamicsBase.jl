module AstrodynamicsBase

    using LinearAlgebra

    include("trigonometry.jl")
    include("body_properties.jl")
    include("elements.jl")
    include("canonical_thrust.jl")

    export acos_safe, asin_safe
    export get_body_soi, get_canonical_param
    export kep2cart, cart2kep, kep2mee, mee2kep, mee2cart, cart2mee, cart2poincare
    export dimensional2canonical_thrust, dimensional2canonical_mdot, dimensional2canonical_isp

end # module
