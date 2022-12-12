"""
Test for Kustaanheimo-Stiefel transformations
"""

#push!(LOAD_PATH,"../src/")
#using AstrodynamicsBase

include("../src/AstrodynamicsBase.jl")
using Printf
using LinearAlgebra
using DifferentialEquations
using Plots
plotly()

## Conversion tests
mu = 1.0
x_cart = AstrodynamicsBase.kep2cart([1.2, 0.1, 0.5, 2.3, 5.3, 0.2], mu)
@printf("x_cart: %1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e\n", x_cart...)

# convert to Kustaanheimo-Stiefel
x_khsf = AstrodynamicsBase.cart2khsf(x_cart)
@printf("x_khsf: %1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e\n", x_khsf...)

# convert back to Cartesian
x_cart_back = AstrodynamicsBase.khsf2cart(x_khsf)
@printf("x_cart: %1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e\n", x_cart_back...)

diff = x_cart_back - x_cart
@printf("diff: %1.6e, %1.6e, %1.6e, %1.6e, %1.6e, %1.6e\n", diff...)

n_test = 10
for i = 1:n_test
    test_cart = rand(6)
    test_diff = test_cart - AstrodynamicsBase.khsf2cart(AstrodynamicsBase.cart2khsf(test_cart))
    println("Test $i : ", norm(test_diff) < 1e-15)
end

## Propagation test
println("Propagation with Kustaanheimo-Stiefel...")
x_cart = AstrodynamicsBase.kep2cart([1.2, 0.1, 0.5, 0.7, 5.3, 0.2], mu)
period = AstrodynamicsBase.get_period(x_cart, mu)
println("Period: $period")
h0 = mu/norm(x_cart[1:3]) - norm(x_cart[4:6])^2/2
u0 = vcat(
    AstrodynamicsBase.cart2khsf(x_cart),
    h0,   # initial negative of energy
    0.0,  # t0
    1.0   # mass at t0
)
tspan = (0.0, 3period)
params = [0.0, 0.0, 0.0, 1e-3, 1e-4]

# ---------------- KH propagation 1 ---------------- #
tstart = time()
prob_khsf = ODEProblem(AstrodynamicsBase.twobody_khsf!,u0,tspan,params);
sol_khsf  = solve(prob_khsf, Tsit5(), reltol=1e-12, abstol=1e-12);
tend = time()
@printf("Elapsed with KH propagation: %4.4f sec\n", tend-tstart)

# convert KH result to cartesian states
s_eval = LinRange(0, sol_khsf.t[end], 1000)
traj_khsf2cart = zeros(6,length(s_eval))
for (i,s) in enumerate(s_eval)
    traj_khsf2cart[:,i] = AstrodynamicsBase.khsf2cart(sol_khsf(s)[1:8])
end

println("Propagation with Kustaanheimo-Stiefel...")
x_cart = AstrodynamicsBase.kep2cart([1.8, 0.3, 0.2, 0.3, 5.3, 0.2], mu)
period = AstrodynamicsBase.get_period(x_cart, mu)
println("Period: $period")
h0 = mu/norm(x_cart[1:3]) - norm(x_cart[4:6])^2/2
u0 = vcat(
    AstrodynamicsBase.cart2khsf(x_cart),
    h0,   # initial negative of energy
    0.0,  # t0
    1.0   # mass at t0
)
params = [0.0, 0.0, 0.0, 1e-3, 1e-4]

# ---------------- KH propagation 2 ---------------- #
tstart = time()
prob_khsf = ODEProblem(AstrodynamicsBase.twobody_khsf!,u0,tspan,params);
sol_khsf2  = solve(prob_khsf, Tsit5(), reltol=1e-12, abstol=1e-12);
tend = time()
@printf("Elapsed with KH propagation: %4.4f sec\n", tend-tstart)

# convert KH result to cartesian states
s_eval = LinRange(0, sol_khsf2.t[end], 1000)
traj_khsf2cart2 = zeros(6,length(s_eval))
for (i,s) in enumerate(s_eval)
    traj_khsf2cart2[:,i] = AstrodynamicsBase.khsf2cart(sol_khsf2(s)[1:8])
end


# ---------------- Plot ---------------- #
println("Plotting...")
ph = plot(frame_style=:box, xlabel="time, TU", ylabel="h, DU^2/TU^2")
plot!(ph, Array(sol_khsf)[10,:], Array(sol_khsf)[9,:], label="Kustaanheimo-Stiefel")
plot!(ph, Array(sol_khsf2)[10,:], Array(sol_khsf2)[9,:], label="traj 2")

ptraj = plot(frame_style=:box, camera=(70,10), xlabel="x, DU", ylabel="y, DU", zlabel="z, DU")
plot!(ptraj, traj_khsf2cart[1,:], traj_khsf2cart[2,:], traj_khsf2cart[3,:], label="Kustaanheimo-Stiefel")
plot!(ptraj, traj_khsf2cart2[1,:], traj_khsf2cart2[2,:], traj_khsf2cart2[3,:], label="traj 2")

cs = palette([:blue, :orange], 4)
pkh_p = plot(frame_style=:box, xlabel="time, TU",)
for i = 1:4
    plot!(pkh_p, sol_khsf.t, Array(sol_khsf)[i,:], label="p_$i", c=cs[i])
    plot!(pkh_p, sol_khsf2.t, Array(sol_khsf2)[i,:], label="p_$i", linestyle=:dash, c=cs[i])
end

pkh_pprime = plot(frame_style=:box, xlabel="time, TU",)
for i = 1:4
    plot!(pkh_pprime, sol_khsf.t, Array(sol_khsf)[4+i,:], label="p'_$i", c=cs[i])
    plot!(pkh_pprime, sol_khsf2.t, Array(sol_khsf2)[4+i,:], label="p'_$i", linestyle=:dash, c=cs[i])
end

display(plot(ptraj, ph, pkh_p, pkh_pprime; layout=(2,2), size=(800,800)))
