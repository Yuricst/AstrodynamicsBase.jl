"""
Test for EOMs
"""

using BenchmarkTools: @btime, @benchmark
push!(LOAD_PATH, "../src")
using AstrodynamicsBase

mu = 398600.44
j2 = 0.00108263
re = 6378.0
params = [mu, j2, re]

oe_kep = [6500, 0.04, deg2rad(93), 0.2, 0.1, 0.5]
oe_mee = AstrodynamicsBase.kep2mee(oe_kep)

du = zeros(length(oe_mee))
println(du)
AstrodynamicsBase.twobody_mee_j2!(du, oe_mee, params, 0.0)
println(du)

println("Benchmarking...")
@benchmark AstrodynamicsBase.twobody_mee_j2!(du, oe_mee, params, 0.0)
