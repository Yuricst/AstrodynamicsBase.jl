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
Ref: See Zimovan thesis 2017 pg.54
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
    inertial2rotating(state_i, theta::Real, om::Real)

Conversion of state-vector from rotating to inertial frame
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
Convert 3-element vector from perifocal to geocentric frame
"""
function perifocal2geocentric(vec::Vector, ω::Real, inc::Real, Ω::Real)
    # rotate by aop
    v1 = _rotmat_ax3(-ω) * vec
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
