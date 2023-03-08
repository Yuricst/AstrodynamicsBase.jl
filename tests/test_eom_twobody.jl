"""
Test for EOMs
"""

using DifferentialEquations
using Plots
# push!(LOAD_PATH, "../src")
# using AstrodynamicsBase
include("../src/AstrodynamicsBase.jl")
gr()

mu = 398600.44
j2 = 0.00108263
re = 6378.0
params = [mu, j2, re]

oe_kep = [42164.0, 0.075, deg2rad(43), 0.2, deg2rad(270), 0.0]
oe_cart = AstrodynamicsBase.kep2cart(oe_kep, mu)

# integrate
tspan = (0, 3*86400)
t_eval = LinRange(tspan[1], tspan[2], 200)
prob = ODEProblem(
    AstrodynamicsBase.twobody_cartesian_j2!, oe_cart, tspan, (mu,j2,re),
    method=Tsit5(), reltol=1e-12, abstol=1e-12,
)
sol = DifferentialEquations.solve(prob)
sv_eval = [sol(t) for t in t_eval]

# convert to rotating frame
theta0 = deg2rad(0)
om = 2π/86164
sv_rot = AstrodynamicsBase.inertial2rotating(sv_eval, t_eval, theta0, om)

# convert back to inertial frame 
sv_inr = AstrodynamicsBase.rotating2inertial(sv_rot, t_eval, theta0, om)

# plot
lscale = 1
ptraj = plot(frame_style=:box, size=(800,600))
plot!(ptraj, Array(sol)[1,:]*lscale, Array(sol)[2,:]*lscale, Array(sol)[3,:]*lscale,)
plot!(ptraj, hcat(sv_rot...)[1,:]*lscale, hcat(sv_rot...)[2,:]*lscale, hcat(sv_rot...)[3,:]*lscale,)
plot!(ptraj, hcat(sv_inr...)[1,:]*lscale, hcat(sv_inr...)[2,:]*lscale, hcat(sv_inr...)[3,:]*lscale,)
display(ptraj)