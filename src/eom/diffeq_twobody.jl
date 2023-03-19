"""
Generic two body differential DifferentialEquations
"""

"""
    twobody_cartesian!(du, u, p, t)

Equations of motion for two-body
"""
function twobody_cartesian!(du, u, p, t)
    # unpack arguemnts
    μ = p[1]
    # compute radius
    r = norm(u[1:3])
    # positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    # velocities
    du[4] = -(μ / r^3) * u[1]
    du[5] = -(μ / r^3) * u[2]
    du[6] = -(μ / r^3) * u[3]
end


"""
    twobody_cartesian_j2!(du, u, p, t)

Equations of motion for two-body with J2.
Parameters are `p = [μ, J2, Re]`.
"""
function twobody_cartesian_j2!(du, u, p, t)
    # unpack arguemnts
    μ, J2, Re = p
    # compute radius
    r = norm(u[1:3])
    Re_r2 = (Re/r)^2
    # positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    # velocities
    du[4] = -(μ / r^3) * u[1] * (1 + 1.5J2*Re_r2*(1 - 5*u[3]^2/r^2))
    du[5] = -(μ / r^3) * u[2] * (1 + 1.5J2*Re_r2*(1 - 5*u[3]^2/r^2))
    du[6] = -(μ / r^3) * u[3] * (1 + 1.5J2*Re_r2*(3 - 5*u[3]^2/r^2))
end


"""
    twobody_cartesian_planar3b!(du, u, p, t)

Equations of motion for two-body with planar third body.
Parameters are `p = [μ, mu_3b, a_3b, n_3b, λ_3b0]`.
"""
function twobody_cartesian_planar3b!(du, u, p, t)
    # unpack arguemnts
    μ, mu_3b, a_3b, n_3b, λ_3b0 = p
    # compute radius
    r = norm(u[1:3])
    # third-body terms
    s = a_3b*[cos(λ_3b0 + n_3b*t), sin(λ_3b0 + n_3b*t), 0.0]
    # compute perturbation via Battin's F(q) function
    d = u[1:3] - s
    q = transpose(u[1:3])*(u[1:3] - 2s)/(transpose(s)*s)
    Fq = q*((3+3q+q^2)/(1 + sqrt(1+q)^3))
    du_3b = -mu_3b/norm(d)^3*(u[1:3] + Fq*s)
    # positions
    du[1] = u[4]
    du[2] = u[5]
    du[3] = u[6]
    # velocities
    du[4] = -(μ / r^3) * u[1] + du_3b[1]
    du[5] = -(μ / r^3) * u[2] + du_3b[2]
    du[6] = -(μ / r^3) * u[3] + du_3b[3]
end
