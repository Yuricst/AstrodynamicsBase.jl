"""
    mod_custom(a, n)

Custom modulo function
Ref: https://stackoverflow.com/questions/1878907/how-can-i-find-the-difference-between-two-angles
"""
function mod_custom(a, n)
    return a - floor(a / n) * n
end


"""
    angle_difference(ϕ_fwd::Real, ϕ_bck::Real)

Compute angle difference for periodic angles between 0 and 2π
"""
function angle_difference(ϕ_fwd::Real, ϕ_bck::Real)
    # modulo based
    dϕ = mod_custom((ϕ_bck - ϕ_fwd + π), 2π) - π
    return dϕ
end


"""
    get_sphere_surface(radius::Real=1, nsph::Int=100, center::Vector=[0,0,0])

Create x,y,z vectors for surface-plotting a sphere.

```julia
xsph,ysph,zsph = get_sphere_surface(6373/lstar);
Plots.surface!(p_xyz, xsph, ysph, zsph ; c=:blues)
```

# Arguments
    - `radius::Real`: radius
    - `nsph::Int`: number of points
    - `center::Vector`: center of sphere

# Returns
    - `Tuple`: vectors of x, y, z coordinates
"""
function get_sphere_surface(radius::Real=1.0, nsph::Int=100, center::Vector=[0,0,0])
    u = range(-π, π; length = nsph)
    v = range(0, π; length = nsph)
    x = radius * cos.(u) * sin.(v)' .+ center[1]
    y = radius * sin.(u) * sin.(v)' .+ center[2]
    z = radius * ones(nsph) * cos.(v)' .+ center[3]
    return x,y,z
end
