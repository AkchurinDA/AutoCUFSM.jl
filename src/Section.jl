# --------------------------------------------------
# MATERIAL
# --------------------------------------------------
struct Material{MaterialPropertyType<:Real}
    ID              ::Integer
    E_xx            ::MaterialPropertyType
    E_yy            ::MaterialPropertyType
    v_xx            ::MaterialPropertyType
    v_yy            ::MaterialPropertyType
    G_xy            ::MaterialPropertyType

    function Material(ID::Integer, E_xx::Real, E_yy::Real, v_xx::Real, v_yy::Real, G_xy::Real)
        material_properties = promote(E_xx, E_yy, v_xx, v_yy, G_xy)

        return new{eltype(material_properties)}(ID, material_properties...)
    end
end

# --------------------------------------------------
# NODE
# --------------------------------------------------
struct Node{CoordinateType<:Real}
    ID              ::Integer
    x               ::CoordinateType
    z               ::CoordinateType
    T               ::Real

    function Node(ID::Integer, x::Real, z::Real, T::Real)
        nodal_coordinates = promote(x, z)

        return new{eltype(nodal_coordinates)}(ID, nodal_coordinates..., T)
    end
end

# --------------------------------------------------
# ELEMENT
# --------------------------------------------------
struct Element{CoordinateType<:Real, MaterialPropertyType<:Real}
    ID              ::Integer
    node_i          ::Node{CoordinateType}
    node_j          ::Node{CoordinateType}
    material        ::Material{MaterialPropertyType}
    t               ::Real
    a
    θ

    function Element(ID::Integer, node_i::Node{T}, node_j::Node{T}, material::Material{S}, t::Real) where {T<:Real, S<:Real}
        x_i, z_i = node_i.x, node_i.z
        x_j, z_j = node_j.x, node_j.z
        dx = x_j - x_i
        dz = z_j - z_i
        a = sqrt(dx ^ 2 + dz ^ 2)
        θ = atan(dz, dx)

        new{T, S}(ID, node_i, node_j, material, t, a, θ)
    end
end

function Element(ID::Integer, node_i::Node{<:Real}, node_j::Node{<:Real}, material::Material, t::Real)
    x_i, z_i = node_i.x, node_i.z
    x_j, z_j = node_j.x, node_j.z

    nodal_coordinates = promote(x_i, z_i, x_j, z_j)

    updated_node_i = Node(node_i.ID, nodal_coordinates[1:2]..., node_i.T)
    updated_node_j = Node(node_j.ID, nodal_coordinates[3:4]..., node_j.T)

    return Element(ID, updated_node_i, updated_node_j, material, t)
end

# --------------------------------------------------
# SECTION
# --------------------------------------------------
struct Section{NN, NE, NM}
    nodes           ::StaticArrays.SVector{NN,     Node}
    elements        ::StaticArrays.SVector{NE,  Element}
    materials       ::StaticArrays.SVector{NM, Material}

    function Section(nodes::StaticArrays.SVector{NN, <:Node}, elements::StaticArrays.SVector{NE, <:Element}, materials::StaticArrays.SVector{NM, <:Material}) where {NN,NE,NM}
        return new{NN,NE,NM}(nodes, elements, materials)
    end
end

function Section(nodes::AbstractVector{<:Node}, elements::AbstractVector{<:Element}, materials::AbstractVector{<:Material})
    NN = length(    nodes)
    NE = length( elements)
    NM = length(materials)

    return Section(
        StaticArrays.SVector{NN}(nodes    ),
        StaticArrays.SVector{NE}(elements ),
        StaticArrays.SVector{NM}(materials))
end