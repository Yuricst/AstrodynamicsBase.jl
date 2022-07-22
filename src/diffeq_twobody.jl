"""
Generic two body differential DifferentialEquations
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