"""
Test for dynamics_mee for rdv with free final time
"""

using LinearAlgebra
import AstrodynamicsBase
using joptimise
using Plots
using DifferentialEquations

include("../src/POLT.jl")

# modify paths as necessary!
spice_dir = ENV["SPICE"]

# get spice kernels
AstrodynamicsBase.spice_furnsh(joinpath(spice_dir, "lsk", "naif0012.tls"))
AstrodynamicsBase.spice_furnsh(joinpath(spice_dir, "spk", "de440.bsp"))

# nondimensional parameters
μ = 1.0
lstar, tstar, vstar = AstrodynamicsBase.get_canonical_param()

# define spacecraft parameters
mstar = 1500   # kg
isp = 3800     # seconds
tmax = 0.33     # Newtons
g0 = 9.81
mdot = tmax/(isp*g0)

# get Earth-Mars positions
et0 = AstrodynamicsBase.spice_utc2et("2024-01-01T00:00:00.000")
tf = [1.2, 1.8]*365*86400/tstar
state0_cart = AstrodynamicsBase.spice_spkssb(
    et0, 3, "ECLIPJ2000", lstar, vstar
)
state_target_cart = AstrodynamicsBase.spice_spkssb(
    et0, 2, "ECLIPJ2000", lstar, vstar
)
state0 = AstrodynamicsBase.cart2mee(state0_cart, μ)
state_target = AstrodynamicsBase.cart2mee(state_target_cart, μ)[1:5]
mass0 = 1.0
println("state0: $state0")
println("state_target: $state_target")

# compute non-dimensional parameters
tmax_nd = tmax * (1/mstar)*(tstar^2/(1e3*lstar))
mdot_nd = mdot *(tstar/mstar)
isp_nd  = isp / tstar

# spacecraft parameters
c1 = tmax_nd
c2 = mdot_nd
ϵ = 1.0

## run minimizer with IPOPT
sn_options = Dict(
    "max_iter" => 100,   # 1500 ~ 2500
    "print_level" => 5,
    "tol" => 1e-3
)

# create problem
use_λ0_hypersphere = false

if use_λ0_hypersphere == true
    x0 = vcat(0.1*rand(8), tf[2])
else
    x0 = vcat(0.1*rand(7), tf[2])
end
println("x0: $x0")

prob = POLT.MEEProblem(
    state0=state0,
    mass0=mass0,
    state_target=state_target,
    tf=tf,
    μ=μ,
    c1=c1,
    c2=c2,
    solver="snopt",
    options=sn_options,
    use_λ0_hypersphere=use_λ0_hypersphere,
    initial_guess=x0,
    method=Tsit5(),
    reltol=1e-11,
    abstol=1e-12,
)
# solve ocp
xopt, fopt, info = POLT.solve_ocp(prob)
tf_opt = xopt[end]

# get cartesian state history
sol_solved = POLT.get_mee_ode_sol(
    xopt, state0, mass0, tf_opt, μ, c1, c2, ϵ,
    false,
    Tsit5(), 1e-11, 1e-12
)
delts, Hf, λ0_check, λLf, λmf = POLT.check_mee_sol(
    sol_solved, state_target, xopt,
    μ, c1, c2, ϵ,
    use_λ0_hypersphere
)
println("Final λL: ", λLf)
println("Final λm: ", λmf)
println("Final Hamiltonian: ", Hf)
if use_λ0_hypersphere
    println("Initial costates norm: ", λ0_check)
end
println("Elements difference: \n", delts)

# -------------------------------------------------------------------------- #
## plots
ts_eval = LinRange(sol_solved.t[1], sol_solved.t[end], 500)
traj_coord = zeros(length(ts_eval), 14)
u_history  = zeros(length(ts_eval))
H_history  = zeros(length(ts_eval))
for (i,t) in enumerate(ts_eval)
    traj_coord[i,:] = vcat(
        AstrodynamicsBase.mee2cart(sol_solved(t)[1:6], μ),
        sol_solved(t)[7:14]
    )
    u_history[i] = POLT.mee_control_magnitude_aposteriori(
        sol_solved(t), 1.0, μ, ϵ, c1, c2
    )
    H_history[i] = POLT.hamiltonian_mee(sol_solved(t)[1:14], [μ, c1, c2, ϵ, 1.0])
end

# plots
xtarget_0_cart = AstrodynamicsBase.mee2cart(vcat(state_target,0.0),μ)
statef_cart = AstrodynamicsBase.keplerder_nostm(
    μ, xtarget_0_cart, 0.0, tf_opt, 1e-12, 25
)
# initial and final orbit
initial_orbit = AstrodynamicsBase.get_orbit(state0, μ, "mee")
final_orbit = AstrodynamicsBase.get_orbit(vcat(state_target,0.0), μ, "mee")

ptraj = plot(
    frame_style=:box, size=(700,500), aspect_ratio=:equal, xlabel="x", ylabel="y"
)
plot!(ptraj, traj_coord[:,1], traj_coord[:,2], linewidth=1.2, label="Trajectory")
scatter!(ptraj, [traj_coord[1,1],], [traj_coord[1,2],], marker=:circle, label="Initial state")
plot!(ptraj, initial_orbit[1,:], initial_orbit[2,:], label="Initial orbit", linestyle=:dash)
plot!(ptraj, final_orbit[1,:], final_orbit[2,:], label="Final orbit", linestyle=:dash)

pcontrol = plot(
    frame_style=:box, size=(700,500), xlabel="Time", ylabel="Control", ylim=[-0.1,1.1]
)
plot!(pcontrol, ts_eval, u_history, linewidth=2.0, label="u")
plot!(pcontrol, ts_eval, traj_coord[:,7], linewidth=2.0, label="m")
plot!(pcontrol, ts_eval, H_history, linewidth=2.0, label="H")
savefig(plot(ptraj, pcontrol), "example_earthmars_mee_freetf.png")
display(plot(ptraj, pcontrol))

println("Done")
