# --------------------------------------------------
# MATERIAL
# --------------------------------------------------
struct Material{MaterialPropertyType<:Real}
    ID          ::Integer
    E_xx        ::MaterialPropertyType
    E_yy        ::MaterialPropertyType
    v_xx        ::MaterialPropertyType
    v_yy        ::MaterialPropertyType
    G_xy        ::MaterialPropertyType

    function Material(ID::Integer, E_xx::Real, E_yy::Real, v_xx::Real, v_yy::Real, G_xy::Real)
        material_properties = promote(E_xx, E_yy, v_xx, v_yy, G_xy)
        return new{eltype(material_properties)}(ID, material_properties...)
    end
end

# --------------------------------------------------
# NODE
# --------------------------------------------------
struct Node{CoordinateType<:Real}
    ID          ::Integer
    x           ::CoordinateType
    z           ::CoordinateType
    T           ::Real

    function Node(ID::Integer, x::Real, z::Real, T::Real)
        nodal_coordinates = promote(x, z)
        return new{eltype(nodal_coordinates)}(ID, nodal_coordinates..., T)
    end
end

# --------------------------------------------------
# ELEMENT
# --------------------------------------------------
struct Element{ElementPropertyType<:Real}
    ID          ::Integer
    node_i      ::Node
    node_j      ::Node
    material    ::Material
    t           ::Real
    # TODO: Define the element properties of interest that should be computed by the inner constructor.
    
    function Element(ID::Integer, node_i::Node{T}, node_j::Node{T}, material::Material, t::Real) where {T <: Real}
        # TODO: Compute geometric properties of the element.
        # TODO: Compute stiffness properties of the element.
       new{typeof(t)}(ID, node_i, node_j, material, t)
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
    # NN = number of nodes
    # NE = number of elements
    # NM = number of materials
    nodes       ::StaticArrays.SVector{NN,     Node}
    elements    ::StaticArrays.SVector{NE,  Element}
    materials   ::StaticArrays.SVector{NM, Material}
end

function Section(nodes::AbstractVector{<:Node}, elements::AbstractVector{<:Element}, materials::AbstractVector{<:Material})
    return Section(
        StaticArrays.SVector{length(    nodes)}(    nodes), 
        StaticArrays.SVector{length( elements)}( elements), 
        StaticArrays.SVector{length(materials)}(materials))
end