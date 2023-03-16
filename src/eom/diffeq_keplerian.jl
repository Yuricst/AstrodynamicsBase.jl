"""
Gauss equation with Keplerian elements
"""


"""
    gauss_keplerian_3b!(du,oe,params,t)

Gauss equation in Keplerian elements 3rd body effect.
Parameters are `p = [mu, mu_3b, a_3b, n_3b, λ_3b0]`.
"""
function gauss_keplerian_3b!(du,oe,params,t)
    # unpack params
    mu, mu_3b, a_3b, n_3b, λ_3b0 = params

    # unpack elements & compute additional parameters
    a,e,i,Ω,ω,f = oe
    cosf = cos(f)
    sinf = sin(f)
    sini = sin(i)
    sinfω = sin(f+ω)
    h = sqrt(a*mu*(1-e^2))
    p = h^2/mu
    rp = a*(1 - e)
    r = p/(1+e*cosf)

    # Gauss perturbation equations
    psi = [
        # multiplies f_r        # multiplies f_theta         # multiplies f_h
        2*a^2/h*e*sinf          2*a^2/h * p/r                0.0
        1/h*p*sinf              1/h*((p+r)*cosf + r*e)       0.0
        0.0                     0.0                          r*cos(f+ω)/h
        0.0                     0.0                          r*sinfω/(h*sini)
       -p*cosf/(e*h)            (p+r)*sinf/(e*h)            -r*sinfω*cos(i)/(h*sini)
        p*cosf/(e*h)            -(p+r)*sinf/(e*h)            0.0
    ]

    # third-body effects, cf. Jagannatha et al 2019
    cosi = cos(i)
    cosθ = cos(ω+f)
    sinθ = sin(ω+f)
    λ_3b = λ_3b0 + n_3b*t
    cos_lmb_Ω = cos(λ_3b - Ω)
    sin_lmb_Ω = sin(λ_3b - Ω)
    rvec = [r,0,0]
    r3b = a_3b*[
        cosθ*cos_lmb_Ω + cosi*sinθ*sin_lmb_Ω
       -sinθ*cos_lmb_Ω + cosi*cosθ*sin_lmb_Ω
       -sini*sin_lmb_Ω
    ]
    rrel = r3b - rvec
    third_body = mu_3b*(rrel/norm(rrel)^3 - r3b/norm(r3b)^3)
    du[1:6] = psi*third_body + [0.0; 0.0; 0.0; 0.0; 0.0; h/r^2]
    return
end
