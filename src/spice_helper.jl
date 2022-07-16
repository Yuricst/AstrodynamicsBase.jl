"""
Spice helpers
"""



"""
    spice_furnsh(loadfle::String)

Wraps furnsh(), furnishes spcie file.
"""
function spice_furnsh(loadfle::String)
    return furnsh(loadfle)
end


"""
    spice_utc2et(str::String)

Wrapper for utc2et()

# Arguments
    - `str::String`: ephemeris, in format "yyyy-mm-ddTHH:MM:SS.fff"
"""
function spice_utc2et(str::String)
    return utc2et(str)
end


"""
    spice_et2utc(et::Float64, format::String='C', prec::Int=3)

Wrapper for et2utc()
"""
function spice_et2utc(et::Float64, format::String="C", prec::Int = 3)
    return et2utc(et, format, prec)
end


"""
	get_spice_locate_function(
		body_id::Int,
		et0::Float64, 
		etf::Float64,
		canonical_params=nothing,
		frame::String="ECLIPJ2000"
	)

Create function with signature:
`locate_spice(epoch) -> r::Vector{Float64}, v::Vector{Float64}`
"""
function get_spice_locate_function(
	body_id::Int,
	et0::Float64, 
	frame::String="ECLIPJ2000"
)
	# if canonical parameters is not provided, create Sun-Earth based ones
	if isnothing(canonical_params)
		# dynamics constants
		MU_SUN = 1.3271244004193938e11
		lstar = 1.495978707e8  # 1AU in km
		vstar = sqrt(MU_SUN / lstar)
		tstar = lstar / vstar
		canonical_params = SunEarthCanonical(
			lstar, vstar, tstar, et0
		)
	end

	# construct functions to locate positions
	function locate_spice(epoch::Float64)
		sv = spkssb(body_id, et0 + epoch*canonical_params.tstar, frame)
		return sv[1:3]/lstar, sv[4:6]/canonical_params.vstar
	end
	return locate_spice, canonical_params
end


"""
	get_spice_locate_function(
		body_id::Int,
		et0_str::String="2022-01-01T00:00:00.00", 
		canonical_params=nothing,
		frame::String="ECLIPJ2000"
	)

Alias with string-based inputs for earliest epoch
"""
function get_spice_locate_function(
	body_id::Int,
	et0_str::String="2022-01-01T00:00:00.00", 
	canonical_params=nothing,
	frame::String="ECLIPJ2000"
)
	# convert epochs
	et0 = utc2et(et0_str)
	return get_spice_locate_function(
		body_id, et0, canonical_params, frame
	)
end



"""
    interp_spkssb(et0, et1, bodyID, frame, lstar, vstar, et_scale=1.0, dt=432000.0)

Create interpolation of body state from SPICE data, between et0 and et1

# Arguments
    - `et0::Float64`: earliest epoch for interpolation
    - `et1::Float`: latest epoch for interpolation
    - `bodyID::Int`: NAIF body ID to interpolate state
    - `frame::String`: SPICE frame to be used, should be "ECLIPJ2000"
    - `lstar::Float`: canonical length scale
    - `vstar::Float`: canonical velocity scale
    - `et_scale::Float`: scaling for time, may be canonical time scale
    - `dt::Float`: time-step to be used, in seconds, recommends around 5 days
    - `verbosity::Int`: level of verbosity, default is 0

# Returns
    - (tuple): Spline1D objects, in order (cxs, cys, czs, cvxs, cvys, cvzs)
"""
function interp_spkssb(
    et0::Float64,
    et1::Float64,
    bodyID::Int,
    frame::String,
    lstar::Real,
    vstar::Real,
    et_scale::Real = 1.0,
    dt::Real = 432000.0,
    verbosity::Int = 0,
)
    # normalize bound epochs
    et0_scaled = et0 / et_scale
    et1_scaled = et1 / et_scale
    dt_scaled = dt / et_scale

    # number of data points to use
    n = Int(ceil((et1_scaled - et0_scaled) / dt_scaled))
    if verbosity > 0
        println("SPICE interpolation of body $bodyID with $n nodes")
    end

    # construct array of time
    ets = convert(Array, LinRange(et0_scaled, et1_scaled, n))

    # arrays to store
    xs, ys, zs, = zeros(n), zeros(n), zeros(n)
    vxs, vys, vzs = zeros(n), zeros(n), zeros(n)

    for idx = 1:length(ets)   #, et in enumerate(ets):
        state = spice_spkssb(ets[idx] * et_scale, bodyID, frame, lstar, vstar)
        xs[idx] = state[1]
        ys[idx] = state[2]
        zs[idx] = state[3]
        vxs[idx] = state[4]
        vys[idx] = state[5]
        vzs[idx] = state[6]
    end
    # fit spline
    cxs = Spline1D(ets, xs; k = 3)
    cys = Spline1D(ets, ys; k = 3)
    czs = Spline1D(ets, zs; k = 3)
    cvxs = Spline1D(ets, vxs; k = 3)
    cvys = Spline1D(ets, vys; k = 3)
    cvzs = Spline1D(ets, vzs; k = 3)
    return (cxs, cys, czs, cvxs, cvys, cvzs)
end


"""
    shift_eclipj2000_to_bodycenter(sv_eclipj2000::Vector, center::Int, et::Float64)

Shift center of state-vector from ECLIPJ2000 SSB to body.
Velocity is also taken as the difference.
"""
function shift_eclipj2000_to_bodycenter(sv_eclipj2000::Vector, center::Int, et::Float64)
    sv_jup = spkssb(center, et, "ECLIPJ2000")
    sv_shifted = sv_eclipj2000 - sv_jup  #vcat(sv_eclipj2000[1:3] - sv_jup[1:3], sv_eclipj2000[4:6])[:]
    return sv_shifted
end


"""
    spice_spkssb(et::Float64, bodyID::Int=3, frame::String="ECLIPJ2000", lstar::Float64=1., vstar::Float64=1.)

Get state-vector of body, wrapper to spice.spkssb()
Must load ephemeris kernel (.tls file) with spice.furnsh() a priori!

Args:
    `et::Float64`: epoch, in ephemeris seconds
    `bodyID::Int`: NAIF body integer
    `frame::String`: string name of frames, default is 'ECLIPJ2000'
    `lstar::Float64`: canonical length-scale
    `vstar::Float64`: canoniacal velocity-scale
"""
function spice_spkssb(
    et::Float64,
    bodyID::Int,
    frame::String = "ECLIPJ2000",
    lstar::Float64 = 1.0,
    vstar::Float64 = 1.0,
    center_body = 10,
)
    sv = spkssb(bodyID, et, frame)   # evaluate state-vector via spkssb
    # check if center is not sun
    if center_body == 10
        return vcat(sv[1:3] / lstar, sv[4:6] / vstar)
    else
        # shift origin
        sv_shifted = shift_eclipj2000_to_bodycenter(sv, center_body, et)
        return vcat(sv_shifted[1:3] / lstar, sv_shifted[4:6] / vstar)
    end
end

