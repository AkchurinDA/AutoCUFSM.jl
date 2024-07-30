function compute_element_properties(node_i::Node{CoordinateType}, node_j::Node{CoordinateType}, t::Real)::NTuple{7} where {CoordinateType <: Real}
    x_i, z_i = node_i.x, node_i.z
    x_j, z_j = node_j.x, node_j.z

    Δx = x_j - x_i
    Δz = z_j - z_i

    w = sqrt(Δx ^ 2 + Δz ^ 2)

    θ = atan(Δz, Δx)

    A = w * t

    x_c = (x_i + x_j) / 2
    z_c = (z_i + z_j) / 2

    return Δx, Δz, w, θ, A, x_c, z_c
end

function compute_section_properties(nodes::StaticArrays.SVector{NN, <:Node}, elements::StaticArrays.SVector{NE, <:Element})::NTuple{6} where {NN, NE}
    Δx_e  = getfield.(elements, :Δx)
    Δz_e  = getfield.(elements, :Δz)
    A_e   = getfield.(elements, :A)
    x_c_e = getfield.(elements, :x_c)
    z_c_e = getfield.(elements, :z_c)

    A_s = sum(A_e)

    x_c_s = sum(x_c_e .* A_e) / A_s
    z_c_s = sum(z_c_e .* A_e) / A_s

    I_xx = sum(((Δz_e .* Δz_e) / 12 + (z_c_e .- z_c_s) .* (z_c_e .- z_c_s)) .* A_e)
    I_zz = sum(((Δx_e .* Δx_e) / 12 + (x_c_e .- x_c_s) .* (x_c_e .- x_c_s)) .* A_e)
    I_xz = sum(((Δz_e .* Δx_e) / 12 + (z_c_e .- z_c_s) .* (x_c_e .- x_c_s)) .* A_e)

    return A_s, x_c_s, z_c_s, I_xx, I_zz, I_xz
end