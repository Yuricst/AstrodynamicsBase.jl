"""
elements
"""



"""
Compute semi-major axis from Cartesian state
"""
function get_semiMajorAxis(state::Array{<:Real,1}, mu::Float64)
	# Initialize Cartesian Polistion and Velocity
    r = state[1:3]
    v = state[4:6]
    a     = 1.0/(2.0/norm(r) - norm(v)^2/mu)     # Semi-major axis
	return a
end


"""
Compute eccentricity from Cartesian state
"""
function get_eccentricity(state::Array{<:Real,1}, mu::Float64)
	# Initialize Cartesian Polistion and Velocity
    r = state[1:3]
    v = state[4:6]
    h = cross(r, v)  # Angular momentum vector
    p = norm(h)*norm(h)/mu                         # Semi-latus rectum
    a = 1.0/(2.0/norm(r) - norm(v)^2/mu)     # Semi-major axis
    # Numerical stability hack for circular and near-circular orbits
    # Ensures that (1-p/a) is always positive
    if isapprox(a, p, atol=1e-9, rtol=1e-8)
        p = a
    end
    e = sqrt(1 - p/a)    # Eccentricity
	return e
end


"""
Compute inclination from Cartesian state
"""
function get_inclination(state::Array{<:Real,1})
    # decompose state to position and velocity vector
    r = state[1:3]
    v = state[4:6]
    # angular momentum
    h = cross(r, v)
    # inclination
    inc = acos(h[3] / norm(h))
    return inc
end


"""
Compute RAAN in radians from Cartesian state
"""
function get_raan(state::Array{<:Real,1})
    # decompose state to position and velocity vector
    r = state[1:3]
    v = state[4:6]
    # angular momentum
    h = cross(r, v)
    # normal direction of xy-plane
    zdir = [0, 0, 1]
    ndir = cross(zdir, h)
    # compute RAAN
    Ω = atan(ndir[2], ndir[1])
    if Ω < 0.0
        return 2π + Ω
    else
        return Ω
    end
end


"""
	get_aop(state::Array{<:Real,1}, mu::Float64)

Compute AOP in radians from Cartesian state
"""
function get_aop(state::Array{<:Real,1}, mu::Float64)
    # decompose state to position and velocity vector
    r = state[1:3]
    v = state[4:6]
    # angular momentum
    h = cross(r, v)
    # normal direction of xy-plane
    zdir = [0, 0, 1]
    ndir = cross(zdir, h)
    # compute eccentricity vector
    ecc = (1 / mu) * cross(v, h) - r / norm(r)
    if (norm(ecc) != 0) && (norm(ndir) != 0)
        ω = acos_safe(dot(ndir, ecc) / (norm(ndir) * norm(ecc)))
        if ecc[3] < 0
            ω = 2π - ω
        end
    else
        ω = 0.0  #get_raan(state)
    end
    return ω
end


"""
Compute true-anomaly in radians from Cartesian state
"""
function get_trueanomaly(state::Array{<:Real,1}, mu::Float64)
    # decompose state to position and velocity vector
    r = state[1:3]
    v = state[4:6]
    # angular momentum
	h_vec = cross(r, v)
    h = norm(h_vec)
    # eccentricity vector
    e_vec = (1/mu) * cross(v, h_vec) - r/norm(r)
    # radial velocity
    vr = dot(v, r) / norm(r)
	θ = atan(h * vr, h^2 / norm(r) - mu)
    # if dot(r, ecc) / (norm(r) * norm(ecc)) <= 1.0
    #     θ = atan(h * vr, h^2 / norm(r) - mu)
    # else
    #     θ = acos(1.0)
    # end
    # wrap around 2π
    if θ < 0
        return 2π + θ
    else
        return θ
    end
end


#export orbit_period
"""
Compute period from Cartesian state
"""
function get_period(state::Array{<:Real,1}, mu::Float64)
    a = get_semiMajorAxis(state, mu)
    if a > 0
        period = 2 * pi * sqrt(a^3 / mu)
    else
        period = Inf
    end
    return period
end


"""
Compute flight-path angle in radians from Cartesian state
FPA is defined between perpendicular to position vector and velocity vector.
"""
function get_fpa(state::Array{<:Real,1}, mu::Float64)
    h = norm(cross(state[1:3], state[4:6]))
    vperp = h / norm(state[1:3])
    vr = mu / h * norm(get_eccentricity(state, mu)) * sin(get_trueanomaly(state, mu))
    gamma = atan(vr / vperp)
    return gamma
end


function anomaly_true_to_mean(θ::Real, e::Real)
	ea = 2*atan(sqrt((1-e)/(1+e)) * tan(θ/2))
	return anomaly_eccentric_to_mean(ea, e)
end


#export anomaly_eccentric_to_mean
"""
Convert eccentric anomaly into mean anomaly.
Arguments:
- `E::Real`: Eccentric anomaly. [rad] or [deg]
- `e::Real`: Eccentricity. [dimensionless]
- `use_degrees:Bool`: If `true` interpret input will be interpreted as being in degrees, and output will be returned in degrees.
Returns:
- `M::Real`: Mean anomaly. [rad] or [deg]
"""
function anomaly_eccentric_to_mean(E::Real, e::Real; use_degrees::Bool=false)
    # Convert degree input
    if use_degrees == true
        E *= pi/180.0
    end

    # Convert eccentric to mean
    M = E - e*sin(E)

    # Convert degree output
    if use_degrees == true
        M *= 180.0/pi
    end

    return M
end

#export anomaly_mean_to_eccentric
"""
Convert mean anomaly into eccentric anomaly.
Arguments:
- `M::Real`: Mean anomaly. [deg] or [deg]
- `e::Real`: Eccentricity. [dimensionless]
- `use_degrees:Bool`: If `true` interpret input will be interpreted as being in degrees, and output will be returned in degrees.
Returns:
- `E::Real`: Eccentric anomaly. [rad] or [deg]
"""
function anomaly_mean_to_eccentric(M::Real, e::Real; use_degrees::Bool=false)
    # Convert degree input
    if use_degrees == true
        M *= pi/180.0
    end

    # Convert mean to eccentric
    max_iter = 15
    epsilson = eps(Float64)*100.0

    # Initialize starting values
    M = mod(M, 2.0*pi)
    if e < 0.8
        E = M
    else
        E = pi
    end

    # Initialize working variable
    f = E - e*sin(E) - M
    i = 0

    # Iterate until convergence
    while abs(f) > epsilson
        f = E - e*sin(E) - M
        E = E - f / (1.0 - e*cos(E))

        # Increase iteration counter
        i += 1
        if i == max_iter
            error("Maximum number of iterations reached before convergence.")
        end
    end

    # Convert degree output
    if use_degrees == true
        E *= 180.0/pi
    end

    return E
end


"""
	keplerian_to_cartesian(x_oe::Array{<:Real,1}, μ::Real)

Conversion from Keplerian to Cartesian.
Input vector in order `[sma, ecc, inc, Ω, ω, θ]`.
"""
function kep2cart(x_oe::Array{<:Real,1}, μ::Real)
	# unpack
	sma, ecc, inc, Ω, ω, θ = x_oe
    # angular momentum
    if ecc <= 1.0
        h = sqrt(sma * μ * (1 - ecc^2))
    else
        h = sqrt(abs(sma) * μ * (ecc^2 - 1))   # hyperbolic case
    end
    # perifocal vector
    x = (h^2 / μ) * (1 / (1 + ecc * cos(θ))) * cos(θ)
    y = (h^2 / μ) * (1 / (1 + ecc * cos(θ))) * sin(θ)
    vx = (μ / h) * (-sin(θ))
    vy = (μ / h) * (ecc + cos(θ))
    rpf = [x, y, 0.0]
    vpf = [vx, vy, 0.0]

    r_rot3 = _perifocal2geocentric(rpf, ω, inc, Ω)
    v_rot3 = _perifocal2geocentric(vpf, ω, inc, Ω)

    # save inertial state
    state_inr = vcat(r_rot3, v_rot3)[:] # cat(r_rot3, v_rot3, dims=(1,1))
    return state_inr
end



#export sCARTtoOSC
"""
	cartesian_to_keplerian(x::Vector, mu::Real, use_degrees::Bool=false)

Convert Cartesian to Keplerian.
Returns elements in order: `[sma, ecc, inc, OMEGA, omega, θ]`
"""
function cart2kep(x::Array{<:Real,1}, mu::Real, use_degrees::Bool=false)

	sma = get_semiMajorAxis(x, mu)
	ecc = get_eccentricity(x, mu)
	inc = get_inclination(x)
	OMEGA = get_raan(x)  #mod(get_raan(state), 2π)
	omega = get_aop(x, mu)  #mod(get_aop(state, μ), 2π)
	θ = get_trueanomaly(x, mu)  # mod(get_trueanomaly(state, μ), 2π)
	x_oe = [sma, ecc, inc, OMEGA, omega, θ]

    # Convert output to degrees if necessary
    if use_degrees == true
        x_oe[3:6] = x_oe[3:6]*180.0/pi
    end
    return x_oe
end


"""
Convert Keplerian to MEE
"""
function kep2mee(oe_kep::Vector, use_ta::Bool = true)
    # unpack
    sma, ecc, inc, Ω, ω, θ = oe_kep
    # convert
    p = sma * (1 - ecc^2)
    f = ecc * cos(ω + Ω)
    g = ecc * sin(ω + Ω)
    h = tan(inc / 2) * cos(Ω)
    k = tan(inc / 2) * sin(Ω)
    if use_ta == true
        L = Ω + ω + θ
    else
        # compute eccentric or hyperbolic anomaly
        if ecc <= 1.0
            EA = 2 * atan(sqrt(1 - ecc) * tan(θ / 2), sqrt(1 + ecc))
        else
            EA = 2 * atan(sqrt(ecc - 1) * tan(θ / 2), sqrt(ecc + 1))
        end
        L = Ω + ω + EA
    end
    return [p, f, g, h, k, L]
end


"""
    modifiedEquinoctial_to_keplerian(mee_elts::Vector)

Convert Keplerian elements to Modified Equinoctial Elements
Input in order `[p, f, g, h, k, L]`, returned in order `[sma, ecc, inc, OMEGA, omega, θ]`
Ref: https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
"""
function mee2kep(mee_elts::Array{<:Real,1})
    # unpack
    p, f, g, h, k, L = mee_elts
    # convert
    sma = p / (1 - f^2 - g^2)
    ecc = sqrt(f^2 + g^2)
    inc = atan(2 * sqrt(h^2 + k^2), 1 - h^2 - k^2)
    ω = atan(g * h - f * k, f * h + g * k)
    Ω = atan(k, h)
    θ = L - atan(g, h)
    return [sma, ecc, inc, Ω, ω, θ]
end



"""
Convert MEE to Cartesian
"""
function mee2cart(oe_mee::Vector, mu::Real)
    # unpack
	p,f,g,h,k,L = oe_mee
	# intermediate values
	s2 = 1 + h^2 + k^2
	alf2 = h^2 - k^2
	w = 1 + f*cos(L) + g*sin(L)
	r = p/w
	# construct position
	x = r/s2 * (cos(L) + alf2*cos(L) + 2*h*k*sin(L))
	y = r/s2 * (sin(L) - alf2*sin(L) + 2*h*k*cos(L))
	z = 2r/s2 * (h*sin(L) - k*cos(L))
	# construct velocity
	vx = -1/s2*sqrt(mu/p)*(sin(L) + alf2*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + alf2*g)
	vy = -1/s2*sqrt(mu/p)*(-cos(L) + alf2*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + alf2*f)
	vz =  2/s2*sqrt(mu/p)*(h*cos(L) + k*sin(L) + f*h + g*k)
	return [x,y,z,vx,vy,vz]
end



"""
Direct conversion from Cartesian to MEE
https://degenerateconic.com/modified-equinoctial-elements.html
"""
function cart2mee(rv::Array{<:Real,1}, mu::Real)
	r = rv[1:3]
	v = rv[4:6]

	rdv      = dot(r,v)
	rhat     = r/norm(r)
	rmag     = norm(r)
	hvec     = cross(r,v)
	hmag     = norm(hvec)
	hhat     = hvec/norm(hvec)
	vhat     = (rmag*v - rdv*rhat)/hmag
	p        = hmag*hmag / mu
	k        = hhat[1]/(1.0 + hhat[3])
	h        = -hhat[2]/(1.0 + hhat[3])
	kk       = k*k
	hh       = h*h
	s2       = 1.0+hh+kk
	tkh      = 2.0*k*h
	ecc      = cross(v,hvec)/mu - rhat
	fhat = [
		1.0-kk+hh,
		tkh,
		-2.0*k
	]
	ghat = [
		tkh,
		1.0+kk-hh,
		2.0*h
	]
	fhat     = fhat/s2
	ghat     = ghat/s2
	f        = dot(ecc,fhat)
	g        = dot(ecc,ghat)
	L        = atan(rhat[2]-vhat[1],rhat[1]+vhat[2])
	return [p,f,g,h,k,L]
end


"""
Convert Cartesian states to Poincare elements
"""
function cart2poincare(rv::Array{<:Real,1}, mu::Real)
	# get MEEs and Keplerian elements
	oe_mee = cartesian_to_modifiedEquinoctial(rv, mu)
	oe_kep = cartesian_to_keplerian(rv, mu)
	# extract
	p,f,g,h,k,_ = oe_mee
	a,e,i,Ω,ω,θ = oe_kep
	# mean anomaly
	ma = anomaly_true_to_mean(θ, e)
	# compute reoccuring values
	sqrt_1_e2 = sqrt(1-e^2)
	# compute poincare elements
	Λ = sqrt(mu*a)
	ξ = g * sqrt(2Λ/(1 + sqrt_1_e2))
	η = f * sqrt(2Λ/(1 + sqrt_1_e2))
	u = k * sqrt(2Λ*sqrt_1_e2/(1 + cos(i)))
	v = h * sqrt(2Λ*sqrt_1_e2/(1 + cos(i)))
	λ = ma + Ω + ω
	return [Λ,ξ,η,u,v,λ]
end
