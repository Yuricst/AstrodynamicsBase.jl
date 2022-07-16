"""
Coordinate transformations
"""

"""
Rotational matrix about axis-1
"""
function _rotmat_ax1(phi::Float64)
    return [1.0 0.0 0.0; 0.0 cos(phi) sin(phi); 0.0 -sin(phi) cos(phi)]
end


"""
Rotational matrix about axis-2
"""
function _rotmat_ax2(phi::Float64)
    return [cos(phi) 0.0 -sin(phi); 0.0 1.0 0.0; sin(phi) 0.0 cos(phi)]
end


"""
Rotational matrix about axis-3
"""
function _rotmat_ax3(phi::Float64)
    return [cos(phi) sin(phi) 0.0; -sin(phi) cos(phi) 0.0; 0.0 0.0 1.0]
end


"""
    rotating2inertial(state_r, theta::Float64, om::Float64)

Conversion of state-vector from rotating to inertial frame
Ref: See Zimovan thesis 2017 pg.54
"""
function rotating2inertial(state_r::Array{<:Real,1}, theta::Float64, om::Float64)
    # construct transformation matrix
    rotmat = [
        cos(theta) -sin(theta) 0.0 0.0 0.0 0.0
        sin(theta) cos(theta) 0.0 0.0 0.0 0.0
        0.0 0.0 1.0 0.0 0.0 0.0
        -om*sin(theta) -om*cos(theta) 0.0 cos(theta) -sin(theta) 0.0
        om*cos(theta) -om*sin(theta) 0.0 sin(theta) cos(theta) 0.0
        0.0 0.0 0.0 0.0 0.0 1.0
    ]
    # transform
    state_i = rotmat * state_r
    return state_i
end


"""
    inertial2rotating(state_r, theta::Float64, om::Float64)

Conversion of state-vector from rotating to inertial frame
Ref: See Zimovan thesis 2017 pg.54
"""
function inertial2rotating(state_i::Array{<:Real,1}, theta::Float64, om::Float64)
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
function perifocal2geocentric(vec::Array{<:Real,1}, ω::Float64, inc::Float64, Ω::Float64)
    # rotate by aop
    v1 = _rotmat_ax3(-ω) * vec
    # rotate by inclination
    v2 = _rotmat_ax1(-inc) * v1
    # rotate by raan
    v_gec = _rotmat_ax3(-Ω) * v2
    return v_gec
end
