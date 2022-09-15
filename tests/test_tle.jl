"""
Tests for using TLE functions
"""

using Plots
using AstrodynamicsBase

μE = 398600.44

sat_dict = AstrodynamicsBase.tle_from_file("gnss.txt")
sat_qzs = AstrodynamicsBase.query("QZS", sat_dict)

println(sat_qzs["QZS-1 (QZSS/PRN 183)"])

ptraj = plot(frame_style=:box, aspect_ratio=:equal)

for key in keys(sat_qzs)
    kep_elts = AstrodynamicsBase.tle2kep(sat_qzs[key])
    prop = AstrodynamicsBase.get_orbit(kep_elts, μE, "keplerian")
    plot!(ptraj, prop[1,:], prop[2,:], prop[3,:], label=key)
end
display(ptraj)
