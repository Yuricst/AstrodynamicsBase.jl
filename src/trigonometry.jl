"""
Safe acos function
"""
function acos_safe(val::Real)
    return acos(max(-1.0, min(1.0, val)))
end


"""
Safe acos function
"""
function asin_safe(val::Real)
    return asin(max(-1.0, min(1.0, val)))
end
