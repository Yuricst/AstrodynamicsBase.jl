"""
    dimensional2canonical_thrust(tmax::Real, mstar::Real, lstar::Real, tstar::Real)

Convert thrust from dimensional to canonical
"""
function dimensional2canonical_thrust(tmax::Real, mstar::Real, lstar::Real, tstar::Real)
    tmax_nd = tmax * (1/mstar)*(tstar^2/(1e3*lstar))
    return tmax_nd
end


"""
    dimensional2canonical_mdot(mdot::Real, mstar::Real, tstar::Real)

Convert thrust from dimensional to canonical
"""
function dimensional2canonical_mdot(mdot::Real, mstar::Real, tstar::Real)
    mdot_nd = mdot *(tstar/mstar)
    return mdot_nd
end


"""
    dimensional2canonical_isp(isp::Real, tstar::Real)

Convert thrust from dimensional to canonical
"""
function dimensional2canonical_isp(isp::Real, tstar::Real)
    isp_nd  = isp / tstar
    return isp_nd
end
