"""
Coordinate transformations
"""

"""
Rotational matrix about axis-1
"""
function _rotmat_ax1(phi::Real)
    return [1.0 0.0 0.0; 0.0 cos(phi) sin(phi); 0.0 -sin(phi) cos(phi)]
end


"""
Rotational matrix about axis-2
"""
function _rotmat_ax2(phi::Real)
    return [cos(phi) 0.0 -sin(phi); 0.0 1.0 0.0; sin(phi) 0.0 cos(phi)]
end


"""
Rotational matrix about axis-3
"""
function _rotmat_ax3(phi::Real)
    return [cos(phi) sin(phi) 0.0; -sin(phi) cos(phi) 0.0; 0.0 0.0 1.0]
end


"""
    rotating2inertial(state_r, theta::Real, om::Real)

Conversion of state-vector from rotating to inertial frame
Ref: See Zimovan thesis 2017 pg.54 eqn. (3.72)
"""
function rotating2inertial(state_r::Vector, theta::Real, om::Real)
    # construct transformation matrix
    rotmat = [
        cos(theta) -sin(theta) 0.0 0.0 0.0 0.0
        sin(theta)  cos(theta) 0.0 0.0 0.0 0.0
        0.0         0.0        1.0 0.0 0.0 0.0
       -om*sin(theta) -om*cos(theta) 0.0 cos(theta) -sin(theta) 0.0
        om*cos(theta) -om*sin(theta) 0.0 sin(theta)  cos(theta) 0.0
        0.0            0.0           0.0 0.0         0.0        1.0
    ]
    # transform
    state_i = rotmat * state_r
    return state_i
end



"""
    rotating2inertial(states_r::Vector, times::Union{Vector,LinRange}, theta0::Real, om::Real)

Conversion of state history from rotating to inertial frame
Ref: See Zimovan thesis 2017 pg.54
"""
function rotating2inertial(states_r::Vector, times::Union{Vector,LinRange}, theta0::Real, om::Real)
    # construct transformation matrix
    states_i = Vector[]
    for (idx,t) in enumerate(times)
        theta = theta0 + om*t
        rotmat = [
            cos(theta) -sin(theta) 0.0 0.0 0.0 0.0
            sin(theta)  cos(theta) 0.0 0.0 0.0 0.0
            0.0         0.0        1.0 0.0 0.0 0.0
           -om*sin(theta) -om*cos(theta) 0.0 cos(theta) -sin(theta) 0.0
            om*cos(theta) -om*sin(theta) 0.0 sin(theta)  cos(theta) 0.0
            0.0            0.0           0.0 0.0         0.0        1.0
        ]
        # transform
        push!(states_i, rotmat * states_r[idx])
    end
    return states_i
end



"""
    inertial2rotating(state_i::Vector, theta::Real, om::Real)

Conversion of state-vector from inertial to rotating frame
Ref: See Zimovan thesis 2017 pg.54
"""
function inertial2rotating(state_i::Vector, theta::Real, om::Real)
    # construct transformation matrix
    rotmat = inv(
        [
            cos(theta) -sin(theta) 0.0 0.0 0.0 0.0
            sin(theta) cos(theta) 0.0 0.0 0.0 0.0
            0.0 0.0 1.0 0.0 0.0 0.0
            -om*sin(theta) -om*cos(theta) 0.0 cos(theta) -sin(theta) 0.0
            om*cos(theta) -om*sin(theta) 0.0 sin(theta) cos(theta) 0.0
            0.0 0.0 0.0 0.0 0.0 1.0
        ],
    )
    # transform
    state_r = rotmat * state_i
    return state_r
end



"""
    inertial2rotating(states_i::Vector, times::Vector, theta0::Real, om::Real)

Conversion of state history from inertial to rotating frame
Ref: See Zimovan thesis 2017 pg.54
"""
function inertial2rotating(states_i::Vector, times::Union{Vector,LinRange}, theta0::Real, om::Real)
    # construct transformation matrix
    states_r = Vector[]
    for (idx,t) in enumerate(times)
        theta = theta0 + om*t
        rotmat = inv(
            [
                cos(theta) -sin(theta) 0.0 0.0 0.0 0.0
                sin(theta) cos(theta) 0.0 0.0 0.0 0.0
                0.0 0.0 1.0 0.0 0.0 0.0
                -om*sin(theta) -om*cos(theta) 0.0 cos(theta) -sin(theta) 0.0
                om*cos(theta) -om*sin(theta) 0.0 sin(theta) cos(theta) 0.0
                0.0 0.0 0.0 0.0 0.0 1.0
            ],
        )
        # transform
        push!(states_r, rotmat * states_i[idx])
    end
    return states_r
end



"""
Convert 3-element vector from perifocal to geocentric frame
"""
function perifocal2geocentric(vec_pf::Vector, ω::Real, inc::Real, Ω::Real)
    # rotate by aop
    v1 = _rotmat_ax3(-ω) * vec_pf
    # rotate by inclination
    v2 = _rotmat_ax1(-inc) * v1
    # rotate by raan
    v_gec = _rotmat_ax3(-Ω) * v2
    return v_gec
end



"""
Convert from Planet-Moon to Sun-Planet CR3BP
"""
function planetmoon2sunplanet(
    state::Vector, μ_moon::Real, μ_sun::Real, theta::Real, om::Real,
    scale_l::Real, scale_v::Real, scale_t::Real
)
    # 1. shift to PLanet-centered, planet-moon rotating frame
    state_1 = state + [μ_moon, 0, 0, 0, 0, 0]
    # 2. convert to PLanet-centered, inertial
    state_2 = rotating2inertial(state_1, theta, om)
    # 3. re-scale
    state_3 = vcat(
        state_2[1:3] * scale_l,
        state_2[4:6] * scale_v
    )
    # 4. convert to Planet-centered, Sun-planet rotating frame
    theta_sp = scale_t * theta
    state_4 = inertial2rotating(state_3, theta_sp, om)
    # 5. shift to Sun-planet barycenter centered, Sun-planet rotating frame
    return state_4 + [1 - μ_sun, 0, 0, 0, 0, 0]
end



# """
#     lvlh2cart(vec_lvlh::Vector, mu::Real, state0::Vector)

# Convert LVLH to Cartesian frame
# """
# function lvlh2cart(vec_lvlh::Vector, mu::Real, state0::Vector)
#     # conversion matrix following Markley and Crassidis 2014 pg.36
#     r, v = state0[1:3], state0[4:6]
#     o3I = -r / norm(r)
#     o2I = -cross(r, v) / norm(cross(r, v))
#     o1I = cross(o2I, o3I)
#     A_IO = reshape(vcat(o1I, o2I, o3I), 3, 3)
#     dv_inertial = A_IO * vec_lvlh
#     return dv_inertial
# end


"""
    eci2lvlh(x_I::Vector, y_B::Vector)

Convert 3-component `y_B` vector in ECI to LVLH using current
position and velocity vector in `x_I`
"""
function eci2lvlh(x_I::Vector, y_B::Vector)
    x_hat = x_I[1:3]/norm(x_I[1:3])
    z_hat = cross(x_I[1:3], x_I[4:6])/norm(cross(x_I[1:3], x_I[4:6]))
    y_hat = cross(z_hat, x_hat)

    R_B2I = [x_hat y_hat z_hat]

    return R_B2I' * y_B
end


"""
    lvlh2eci(x_I::Vector, y_B::Vector)
    
Convert 3-component `y_B` vector in LVLH to ECI using current
position and velocity vector in `x_I`
"""
function lvlh2eci(x_I::Vector, y_B::Vector)
    x_hat = x_I[1:3]/norm(x_I[1:3])
    z_hat = cross(x_I[1:3], x_I[4:6])/norm(cross(x_I[1:3], x_I[4:6]))
    y_hat = cross(z_hat, x_hat)

    R_B2I = [x_hat y_hat z_hat]

    return R_B2I * y_B
end


"""
Convert angle in arc second to radians
"""
function arcsec2rad(arcsec)
    return arcsec * pi / (180 * 3600);
end


"""
Convert state in ECI to ECEF

# Arguments 
    - `jd`: Julian date, e.g. 2459449.5
    - `x_eci`: State vector in ECI coordinates, may be just position (length 3) or position & velocity (length 6)

# Returns 
    State vector in ECEF coordinates
"""
function eci2ecef(jd::Real, x_eci::Vector, omega_e::Real=7.29211585275553e-005)
    # Compute the Julian centuries from J2000.0 epoch
    T = (jd - 2451545.0) / 36525.0

    # Compute the mean obliquity of the ecliptic in radians
    eps0 = deg2rad(23.439291 - 0.0130042 * T - 1.64e-7 * T^2 + 5.04e-7 * T^3)

    # Compute the nutation in obliquity in radians
    d_eps = arcsec2rad(9.20 * sind(0.0 * T) + 0.57 * sind(2.0 * T) + 0.10 * sind(2.0 * T - 2.0 * T) - 0.09 * sind(0.0 * T)) / 3600.0

    # Compute the nutation in longitude in radians
    d_psi = arcsec2rad(0.00 * cosd(0.0 * T) + 0.16 * cosd(2.0 * T) + 0.00 * cosd(2.0 * T - 2.0 * T) + 0.00 * cosd(0.0 * T)) / 3600.0

    # Compute the true obliquity of the ecliptic in radians
    eps = eps0 + d_eps

    # Compute the Greenwich Mean Sidereal Time (GMST) in radians
    gmst_rad = deg2rad(280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T^2 - T^3 / 38710000.0 + d_psi * cos(eps))

    # Define the rotation matrix from ECI to ECEF coordinates
    R_eci_to_ecef = [cos(gmst_rad) sin(gmst_rad) 0; -sin(gmst_rad) cos(gmst_rad) 0; 0 0 1]

    # Compute the position vector in ECEF coordinates
    r_ecef = R_eci_to_ecef * x_eci[1:3]

    if length(x_eci) == 3
        return r_ecef
    else
        # Compute the velocity vector in ECEF coordinates
        omega_e_vec = [0 0 omega_e]  # Earth's rotation rate vector in ECEF coordinates
        v_ecef = R_eci_to_ecef * v_eci - cross(omega_e_vec, r_ecef)
        return vcat(r_ecef, v_ecef)
    end
end
