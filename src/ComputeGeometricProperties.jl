function compute_element_properties(node_i::Node{CoordinateType}, node_j::Node{CoordinateType}, t::Real)::NTuple{3}
    x_i, z_i = node_i.x, node_i.z
    x_j, z_j = node_j.x, node_j.z

    dx = x_j - x_i
    dz = z_j - z_i

    a = sqrt(dx ^ 2 + dz ^ 2)
    A = a * t

    θ = atan(dz, dx)

    return a, A, θ
end

function compute_section_properties()

end
