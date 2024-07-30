# --------------------------------------------------
# MATERIAL
# --------------------------------------------------
struct Material{MaterialPropertyType<:Real}
    E_xx            ::MaterialPropertyType
    E_yy            ::MaterialPropertyType
    v_xx            ::MaterialPropertyType
    v_yy            ::MaterialPropertyType
    G_xy            ::MaterialPropertyType

    function Material(E_xx::Real, E_yy::Real, v_xx::Real, v_yy::Real, G_xy::Real)
        material_properties = promote(E_xx, E_yy, v_xx, v_yy, G_xy)

        return new{eltype(material_properties)}(material_properties...)
    end
end

# --------------------------------------------------
# NODE
# --------------------------------------------------
struct Node{CoordinateType<:Real, StressType<:Real}
    x               ::CoordinateType
    z               ::CoordinateType
    σ               ::StressType

    function Node(x::Real, z::Real, σ::Real)
        nodal_coordinates = promote(x, z)

        return new{eltype(nodal_coordinates), typeof(σ)}(nodal_coordinates..., σ)
    end
end

# --------------------------------------------------
# ELEMENT
# --------------------------------------------------
struct Element{CoordinateType<:Real, MaterialPropertyType<:Real, ElementPropertyType<:Real}
    node_i          ::Node{CoordinateType}
    node_j          ::Node{CoordinateType}
    material        ::Material{MaterialPropertyType}
    t               ::ElementPropertyType
    Δx              ::ElementPropertyType
    Δz              ::ElementPropertyType
    w               ::ElementPropertyType
    θ               ::ElementPropertyType
    A               ::ElementPropertyType
    x_c             ::ElementPropertyType
    z_c             ::ElementPropertyType

    function Element(node_i::Node{CoordinateType}, node_j::Node{CoordinateType}, material::Material{MaterialPropertyType}, t::Real) where {CoordinateType<:Real, MaterialPropertyType<:Real}
        element_properties = compute_element_properties(node_i, node_j, t)
        element_properties = promote(t, element_properties...)

        new{CoordinateType, MaterialPropertyType, eltype(element_properties)}(node_i, node_j, material, element_properties...)
    end
end

function Element(node_i::Node{<:Real}, node_j::Node{<:Real}, material::Material, t::Real)
    x_i, z_i = node_i.x, node_i.z
    x_j, z_j = node_j.x, node_j.z

    nodal_coordinates = promote(x_i, z_i, x_j, z_j)

    updated_node_i = Node(nodal_coordinates[1:2]..., node_i.σ)
    updated_node_j = Node(nodal_coordinates[3:4]..., node_j.σ)

    return Element(updated_node_i, updated_node_j, material, t)
end

# --------------------------------------------------
# SECTION
# --------------------------------------------------
struct Section{NN, NE, NM, SectionPropertyType<:Real}
    materials       ::StaticArrays.SVector{NM, Material}
    nodes           ::StaticArrays.SVector{NN, Node}
    elements        ::StaticArrays.SVector{NE, Element}
    A               ::SectionPropertyType
    x_c             ::SectionPropertyType
    z_c             ::SectionPropertyType
    I_xx            ::SectionPropertyType
    I_zz            ::SectionPropertyType
    I_xz            ::SectionPropertyType

    function Section(materials::StaticArrays.SVector{NM, <:Material}, nodes::StaticArrays.SVector{NN, <:Node}, elements::StaticArrays.SVector{NE, <:Element}) where {NM, NN, NE}
        section_properties = compute_section_properties(nodes, elements)
        section_properties = promote(section_properties...)

        return new{NN, NE, NM, eltype(section_properties)}(materials, nodes, elements, section_properties...)
    end
end

Section(materials::AbstractVector{<:Material}, nodes::AbstractVector{<:Node}, elements::AbstractVector{<:Element}) = Section(StaticArrays.SVector{length(materials)}(materials), StaticArrays.SVector{length(nodes)}(nodes), StaticArrays.SVector{length(elements)}(elements))