"""
Functions for using Kustaanheimo-Stiefel transformations
"""


"""
    cart2khsf(x_cart::Array{<:Real,1}, p1::Real=0.1)

Convert Cartesian state to Kustaanheimo-Sfiefel state
See:
    Tracy and Manchester. Low-Thrust Trajectory Optimization using
    the Kustaanheimo-Stiefel Transformation. 31st AIAA/AAS Space Flight
    Mechanics Meeting in Charlotte, North Carolina, 2021.
"""
function cart2khsf(x_cart::Array{<:Real,1}, p1::Real=0.1)
    # convert position
    absx_x1 = norm(x_cart[1:3]) + x_cart[1]
    p4 = sqrt(0.5*absx_x1 - p1^2)
    p2 = (x_cart[2]*p1 + x_cart[3]*p4)/absx_x1
    p3 = (x_cart[3]*p1 - x_cart[2]*p4)/absx_x1
    # convert velocity
    p1_prime = 0.5*( p1*x_cart[4] + p2*x_cart[5] + p3*x_cart[6])
    p2_prime = 0.5*(-p2*x_cart[4] + p1*x_cart[5] + p4*x_cart[6])
    p3_prime = 0.5*(-p3*x_cart[4] - p4*x_cart[5] + p1*x_cart[6])
    p4_prime = 0.5*( p4*x_cart[4] - p3*x_cart[5] + p2*x_cart[6])
    return [
        p1, p2, p3, p4,
        p1_prime, p2_prime, p3_prime, p4_prime
    ]
end


"""
    khsf2cart(p_khsf::Array{<:Real,1})

Convert Kustaanheimo-Sfiefel state to Cartesian state
See:
    Tracy and Manchester. Low-Thrust Trajectory Optimization using
    the Kustaanheimo-Stiefel Transformation. 31st AIAA/AAS Space Flight
    Mechanics Meeting in Charlotte, North Carolina, 2021.
"""
function khsf2cart(p_khsf::Array{<:Real,1})
    # convert position
    x1 = p_khsf[1]^2 - p_khsf[2]^2 - p_khsf[3]^2 + p_khsf[4]^2
    x2 = 2 * (p_khsf[1]*p_khsf[2] - p_khsf[3]*p_khsf[4])
    x3 = 2 * (p_khsf[1]*p_khsf[3] + p_khsf[2]*p_khsf[4])
    # convert velocity
    factor = 2 / (p_khsf[1]^2 + p_khsf[2]^2 + p_khsf[3]^2 + p_khsf[4]^2)
    x1_dot = factor*(
        p_khsf[1]*p_khsf[5] - p_khsf[2]*p_khsf[6] - p_khsf[3]*p_khsf[7] + p_khsf[4]*p_khsf[8]
    )
    x2_dot = factor*(
        p_khsf[2]*p_khsf[5] + p_khsf[1]*p_khsf[6] - p_khsf[4]*p_khsf[7] - p_khsf[3]*p_khsf[8]
    )
    x3_dot = factor*(
        p_khsf[3]*p_khsf[5] + p_khsf[4]*p_khsf[6] + p_khsf[1]*p_khsf[7] + p_khsf[2]*p_khsf[8]
    )
    return [x1, x2, x3, x1_dot, x2_dot, x3_dot]
end



"""
Equations of motion for two-body with perturbation in Kustaanheimo-Stiefel state.
State `u = [p, p', h, t, mass]` where t is dimensional-time, h is negative energy.
Parameters are `params = [τ, θ, β, tmax, mdot]`.
"""
function twobody_khsf!(du, u, params, s)
    # derivative of p
    du[1:4] = u[5:8]
    # accel = τ × tmax / mass * [cθ*cβ, sθ*cβ, sβ, 0]
    accel = params[1] * params[4]/u[11] * [
        cos(params[2]) * cos(params[3]);
        sin(params[2]) * cos(params[3]);
        sin(params[3]);
        0;
    ]
    # accel = params[4]/u[11] * [
    #     params[1],
    #     params[2],
    #     params[3]
    #     0;
    # ]
    # L = [
    #     u[1] -u[2] -u[3]  u[4];
    #     u[2]  u[1] -u[4] -u[3];
    #     u[3]  u[4]  u[1]  u[2];
    #     u[4] -u[3]  u[2] -u[1]
    # ]
    # derivative of p'
    L_T = [
         u[1]  u[2]  u[3]  u[4];
        -u[2]  u[1]  u[4] -u[3];
        -u[3] -u[4]  u[1]  u[2];
         u[4] -u[3]  u[2] -u[1];
    ]
    pnorm2  = norm(u[1:4])^2
    du[5:8] = -u[9]/2*u[1:4] + pnorm2/2 * L_T * accel
    # derivative of h
    du[9]   = - 2 * dot([u[5] u[6] u[7] u[8]], L_T, accel)
    # derivative of t
    du[10]  = pnorm2
    # derivative of mass
    du[11]  = -params[1] * abs(params[5])
end
