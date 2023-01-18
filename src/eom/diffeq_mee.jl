"""
Equations of motion in MEE
Parameters are `params = [μ, τ, θ, β, tmax, mdot]`
"""
function twobody_mee!(du, u, params, t)
    μ, τ, θ, β, tmax, mdot = params
    p,f,g,h,k,L = u
    cosL = cos(L)
    sinL = sin(L)
    w  = 1 + f*cosL + g*sinL
    s2 = 1 + h^2 + k^2
    hsinL_kcosL = (h*sinL - k*cosL)
    sqrt_p_μ = sqrt(p/μ)
    # compute B-matrix
    B = sqrt_p_μ*[
         0    2p/w              0;
         sinL ((1+w)*cosL+f)/w -g/w*hsinL_kcosL;
        -cosL ((1+w)*sinL+g)/w  f/w*hsinL_kcosL;
         0    0                 1/w*s2/2*cosL;
         0    0                 1/w*s2/2*sinL;
         0    0                 1/w*hsinL_kcosL;
    ]
    # copmute D-matrix
    D = [0, 0, 0, 0, 0, sqrt(μ/p^3)*w^2]
    # compute derivatives of MEE and mass
    i_thrust = [
        cos(params[3]) * cos(params[4]);
        sin(params[3]) * cos(params[4]);
        sin(params[4]);
    ]
    du[1:6]  = tmax*τ/u[7] * B*i_thrust + D
    du[7]    = -τ * abs(mdot)
end


"""
Equations of motion in MEE
Parameters are `params = [μ, j2, re]`
"""
function twobody_mee_j2!(du, u, params, t)
    μ, j2, re = params
    p,f,g,h,k,L = u
    cosL = cos(L)
    sinL = sin(L)
    w  = 1 + f*cosL + g*sinL
    s2 = 1 + h^2 + k^2
    hsinL_kcosL = (h*sinL - k*cosL)
    sqrt_p_μ = sqrt(p/μ)
    # compute B-matrix
    B = sqrt_p_μ*[
         0    2p/w              0;
         sinL ((1+w)*cosL+f)/w -g/w*hsinL_kcosL;
        -cosL ((1+w)*sinL+g)/w  f/w*hsinL_kcosL;
         0    0                 1/w*s2/2*cosL;
         0    0                 1/w*s2/2*sinL;
         0    0                 1/w*hsinL_kcosL;
    ]
    # copmute D-matrix
    D = [0, 0, 0, 0, 0, sqrt(μ/p^3)*w^2]
    # compute derivatives of MEE and mass
    r = p/w
    coef_j2 = μ*j2*re^2/r^4
    denom = s2^2
    Δ_J2 = [
        -3/2 * coef_j2 * (1 - 12(h*sinL - k*cosL)^2/denom),
        -12  * coef_j2 * (h*sinL - k*cosL)*(h*cosL + k*sinL)/denom,
        -6   * coef_j2 * (1-h^2-k^2)*(h*sinL - k*cosL)/denom,
    ]
    du .= B*Δ_J2 + D
end
