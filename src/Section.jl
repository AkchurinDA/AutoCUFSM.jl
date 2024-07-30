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
struct Element{CoordinateType<:Real, StressType<:Real, MaterialPropertyType<:Real}
    node_i          ::Node{CoordinateType, StressType}
    node_j          ::Node{CoordinateType, StressType}
    node_i_ID       ::Integer
    node_j_ID       ::Integer
    material        ::Material{MaterialPropertyType}
    material_ID     ::Integer
    t    
    Δx   
    Δz   
    w    
    θ    
    A    
    x_c  
    z_c  

    function Element(
        node_i::Node{CoordinateType, StressType}, node_j::Node{CoordinateType, StressType}, node_i_ID::Integer, node_j_ID::Integer,
        material::Material{MaterialPropertyType}, material_ID::Integer,
        t::Real) where {CoordinateType<:Real, StressType<:Real, MaterialPropertyType<:Real}
        element_properties = compute_element_properties(node_i, node_j, t)
        element_properties = promote(t, element_properties...)

        new{CoordinateType, StressType, MaterialPropertyType}(node_i, node_j, node_i_ID, node_j_ID, material, material_ID, element_properties...)
    end
end

function Element(
    nodes, node_i_ID::Integer, node_j_ID::Integer, 
    materials, material_ID::Integer, 
    t::Real)
    node_i = nodes[node_i_ID]
    node_j = nodes[node_j_ID]

    x_i, z_i = node_i.x, node_i.z
    x_j, z_j = node_j.x, node_j.z

    nodal_coordinates = promote(x_i, z_i, x_j, z_j)

    σ_i = node_i.σ
    σ_j = node_j.σ

    stresses = promote(σ_i, σ_j)

    updated_node_i = Node(nodal_coordinates[1:2]..., stresses[1])
    updated_node_j = Node(nodal_coordinates[3:4]..., stresses[2])

    material = materials[material_ID]

    return Element(updated_node_i, updated_node_j, node_i_ID, node_j_ID, material, material_ID, t)
end

# --------------------------------------------------
# SECTION
# --------------------------------------------------
struct Section{NN, NE, NM}
    materials       ::StaticArrays.SVector{NM, Material}
    nodes           ::StaticArrays.SVector{NN, Node}
    elements        ::StaticArrays.SVector{NE, Element}
    A               
    x_c             
    z_c             
    I_xx            
    I_zz            
    I_xz            

    function Section(materials::StaticArrays.SVector{NM, <:Material}, nodes::StaticArrays.SVector{NN, <:Node}, elements::StaticArrays.SVector{NE, <:Element}) where {NM, NN, NE}
        section_properties = compute_section_properties(nodes, elements)
        section_properties = promote(section_properties...)

        return new{NN, NE, NM}(materials, nodes, elements, section_properties...)
    end
end

Section(materials::AbstractVector{<:Material}, nodes::AbstractVector{<:Node}, elements::AbstractVector{<:Element}) = Section(StaticArrays.SVector{length(materials)}(materials), StaticArrays.SVector{length(nodes)}(nodes), StaticArrays.SVector{length(elements)}(elements))