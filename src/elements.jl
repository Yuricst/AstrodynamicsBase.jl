"""
elements
"""


"""
Convert Keplerian to MEE
"""
function kep2mee(oe_kep::Vector)
    # unpack
    a,e,i,raan,om,ta = oe_kep
    # compute MEEs
    p = a*(1-e^2)
    f = e*cos(raan+om)
    g = e*sin(raan+om)
    h = tan(i/2)*cos(raan)
    k = tan(i/2)*sin(raan)
    l = raan + om + ta
    return [p,f,g,h,k,l]
end


"""
Convert MEE to Cartesian
"""
function mee2cart(oe_mee::Vector, μ::Real)
    # unpack
    p,f,g,h,k,l = oe_mee
    sinL = sin(l)
    cosL = cos(l)
    α2 = h^2 - k^2
    s2 = 1 + h^2 + k^2
    w = 1 + f*cosL + g*sinL
    r = p/w
    sqrt_μp = sqrt(μ/p)
    # state vector in Cartesian
    cart = [
        r/s2  * (cosL + α2*cosL + 2*h*k*sinL),
        r/s2  * (sinL - α2*sinL + 2*h*k*cosL),
        2r/s2 * (h*sinL - k*cosL),
        -1/s2 * sqrt_μp * ( sinL + α2*sinL - 2*h*k*cosL + g - 2*f*h*k + α2*g),
        -1/s2 * sqrt_μp * (-cosL + α2*cosL + 2*h*k*sinL - f + 2*g*h*k + α2*f),
         2/s2 * sqrt_μp * (h*cosL + k*sinL + f*h + g*k),
    ]
    return cart
end