"""
    get_orbit(state, μ::Real, type::String="mee", steps::Int=100)

Get cartesian history of single revolution of orbit.
The returned object is 6-by-n, where n is the time-steps.
"""
function get_orbit(state, μ::Real, type::String="mee", steps::Int=100)
    if type == "mee"
        x0 = mee2cart(state, μ)
    elseif type == "keplerian"
        x0 = kep2cart(state, μ)
    elseif type == "cartesian"
        x0 = state
    end
    period = get_period(x0, μ)
    ts = LinRange(0.0, period, steps)
    _, prop = keplerder_nostm(μ, x0, 0.0, ts)
    return prop
end
