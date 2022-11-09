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
tspan = (0.0, 300.0)
params = [0.0, 0.0, 0.0, 1e-3, 1e-4]

# ---------------- KH propagation ---------------- #
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

# ---------------- MEE propagation ---------------- #
x_mee = AstrodynamicsBase.cart2mee(x_cart, mu)
tstart = time()
prob_mee = ODEProblem(
    AstrodynamicsBase.twobody_mee!,
    vcat(x_mee,1.0),
    [0,sol_khsf.u[end][10]],
    vcat(mu,params),
);
sol_mee  = solve(prob_mee, Tsit5(), reltol=1e-12, abstol=1e-12);
tend = time()
@printf("Elapsed with MEE propagation: %4.4f sec\n", tend-tstart)

# calculate (negative of) specific energy at each MEE state
h_mee = zeros(length(sol_mee.t))
for (i,u) in enumerate(sol_mee.u)
    sv_cart = AstrodynamicsBase.mee2cart(u, mu)
    h_mee[i] = mu/norm(sv_cart[1:3]) - norm(sv_cart[4:6])^2/2
end

# convert KH result to cartesian states
t_eval = LinRange(0, sol_mee.t[end], 1000)
traj_mee2cart = zeros(6,length(t_eval))
for (i,t) in enumerate(t_eval)
    traj_mee2cart[:,i] = AstrodynamicsBase.mee2cart(sol_mee(t)[1:6], mu)
end

# ---------------- Cartesian propagation for comparison ---------------- #
tstart = time()
prob_cart = ODEProblem(AstrodynamicsBase.twobody_cartesian!,x_cart,[0,sol_khsf.u[end][10]],[mu,]);
sol_cart  = solve(prob_cart, Tsit5(), reltol=1e-12, abstol=1e-12);
tend = time()
@printf("Elapsed with Cartesian propagation: %4.4f sec\n", tend-tstart)

# calculate (negative of) specific energy at each Cartesian state
h_cart = zeros(length(sol_cart.t))
for (i,u) in enumerate(sol_cart.u)
    h_cart[i] = mu/norm(u[1:3]) - norm(u[4:6])^2/2
end

# ---------------- Plot ---------------- #
println("Plotting...")
ph = plot(frame_style=:box, xlabel="time, TU", ylabel="h, DU^2/TU^2")
plot!(ph, Array(sol_khsf)[10,:], Array(sol_khsf)[9,:], label="Kustaanheimo-Stiefel")
plot!(ph, sol_cart.t, h_cart, label="Cartesian")
plot!(ph, sol_mee.t, h_mee, label="MEE")

ptraj = plot(frame_style=:box, camera=(70,10), xlabel="x, DU", ylabel="y, DU", zlabel="z, DU")
plot!(ptraj, traj_khsf2cart[1,:], traj_khsf2cart[2,:], traj_khsf2cart[3,:], label="Kustaanheimo-Stiefel")
plot!(ptraj, traj_mee2cart[1,:], traj_mee2cart[2,:], traj_mee2cart[3,:], label="MEE")

display(plot(ptraj, ph; layout=(1,2), size=(800,400)))

