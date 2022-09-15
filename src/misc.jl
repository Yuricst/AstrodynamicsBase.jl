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
