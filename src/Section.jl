struct Section
end

struct Node{CoordinateType<:Real}
    ID::Integer
    x::CoordinateType
    z::CoordinateType
end

Node(ID::Integer, x::Real, z::Real) = Node(ID, promote(x, z)...)

struct Element
    ID::Integer
    NodeI::Node
    NodeJ::Node
end

struct Material{MaterialPropertyType<:Real}
    ID::Integer
    E_xx::MaterialPropertyType
    E_yy::MaterialPropertyType
    v_xx::MaterialPropertyType
    v_yy::MaterialPropertyType
    G_xy::MaterialPropertyType
end

Material(ID::Integer, E_xx::Real, E_yy::Real, v_xx::Real, v_yy::Real, G_xy::Real) = Material(ID, promote(E_xx, E_yy, v_xx, v_yy, G_xy)...)