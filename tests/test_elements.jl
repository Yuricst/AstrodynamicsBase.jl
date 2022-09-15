"""
Test element conversions
"""

using Plots
using LinearAlgebra
using AstrodynamicsBase

mu = 1.0
n_random_test = 100

# RANDOMIZED CASES
println("----------------\nCartesian <-> Keplerian (random)")
failed_cart_kep = []
for i = 1:n_random_test
    x_cart = 2rand(6) - ones(6)
    x_kep = AstrodynamicsBase.cart2kep(x_cart, mu)
    x_cart_back = AstrodynamicsBase.kep2cart(x_kep, mu)
    if (x_cart_back ≈ x_cart) == false
        push!(failed_cart_kep, [x_cart_back, x_cart])
    end
end
if length(failed_cart_kep) == 0
    println("Success!")
else
    println("Failed ", length(failed_cart_kep), " cases!")
end

println("----------------\nCartesian <-> MEE (random)")
failed_cart_mee = []
for i = 1:n_random_test
    x_cart = 2rand(6) - ones(6)
    x_mee = AstrodynamicsBase.cart2mee(x_cart, mu)
    x_cart_back = AstrodynamicsBase.mee2cart(x_mee, mu)
    if (x_cart_back ≈ x_cart) == false
        push!(failed_cart_mee, [x_cart_back, x_cart])
    end
end
if length(failed_cart_mee) == 0
    println("Success!")
end

# EDGE CASES
println("----------------\nEdge 1. Planar prograde")
x_cart = [1.0, 0.1, 0.0, 0.0, 1.2, 0.0]
x_kep = AstrodynamicsBase.cart2kep(x_cart, mu)
x_mee = AstrodynamicsBase.cart2mee(x_cart, mu)
x_kep2cart = AstrodynamicsBase.kep2cart(x_kep, mu)
x_mee2cart = AstrodynamicsBase.mee2cart(x_mee, mu)
if (x_kep2cart ≈ x_cart) == false
    println("Cartesian <-> Keplerian failed! x_kep: $x_kep")
    println(x_cart)
    println(x_kep2cart)
else
    println("Cartesian <-> Keplerian success!")
end
if (x_mee2cart ≈ x_cart) == false
    println("Cartesian <-> MEE failed!")
else
    println("Cartesian <-> MEE success!")
end


println("----------------\nEdge 2. Circular")
x_cart = [1/sqrt(2), 0.0, 1/sqrt(2), 0.0, 1.0, 0.0]
# prop = AstrodynamicsBase.get_orbit(x_cart, 1.0, "cartesian")
# ptraj = plot(prop[1,:], prop[2,:], prop[3,:])
# plot!(ptraj, [x_cart[1]], [x_cart[2]], [x_cart[3]], marker=:circle)
x_kep = AstrodynamicsBase.cart2kep(x_cart, mu)
x_mee = AstrodynamicsBase.cart2mee(x_cart, mu)
x_kep2cart = AstrodynamicsBase.kep2cart(x_kep, mu)
x_mee2cart = AstrodynamicsBase.mee2cart(x_mee, mu)
if (x_kep2cart ≈ x_cart) == false
    println("Cartesian <-> Keplerian failed! x_kep: $x_kep")
    println(x_cart)
    println(x_kep2cart)
else
    println("Cartesian <-> Keplerian success!")
end
if (x_mee2cart ≈ x_cart) == false
    println("Cartesian <-> MEE failed!")
else
    println("Cartesian <-> MEE success!")
end


println("----------------\nEdge 3. Planar retrograde")
x_cart = [1.0, 0.1, 0.0, 0.0, -1.2, 0.0]
prop = AstrodynamicsBase.get_orbit(x_cart, 1.0, "cartesian")
ptraj = plot(prop[1,:], prop[2,:])
plot!(ptraj, [x_cart[1]], [x_cart[2]], marker=:circle)
#vcat(2rand(2)-ones(2), 0.0, 2rand(2)-ones(2), 0.0)
#[0.45720038370837734, -0.07234773411156037, 0.0, 0.4171354636420461, 0.6572732479628429, 0.0]

x_kep = AstrodynamicsBase.cart2kep(x_cart, mu)
x_mee = AstrodynamicsBase.cart2mee(x_cart, mu)
x_kep2cart = AstrodynamicsBase.kep2cart(x_kep, mu)
x_mee2cart = AstrodynamicsBase.mee2cart(x_mee, mu)
if (x_kep2cart ≈ x_cart) == false
    println("Cartesian <-> Keplerian failed! x_kep: $x_kep")
    println(x_cart)
    println(x_kep2cart)
end
if (x_mee2cart ≈ x_cart) == false
    println("Cartesian <-> MEE failed!")
end


display(ptraj)
