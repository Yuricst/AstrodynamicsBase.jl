"""
Test for EOMs
"""

using BenchmarkTools: @btime, @benchmark
push!(LOAD_PATH, "../src")
using AstrodynamicsBase


jd = 2.459604625000000e+06
r_eci = [
    -2.919503698677868e+05,
    -2.213642751580328e+05,
    -6.738911793407053e+04,
]

r_ecef = AstrodynamicsBase.eci2ecef(jd, r_eci)
println("r_ecef = \n$r_ecef")

test_ecef = [
    2.462897069637153e+05,
    2.712536470278437e+05,
    -6.738911793407053e+04,
]
@show r_ecef - test_ecef;

