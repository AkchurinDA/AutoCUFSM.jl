using AutoCUFSM
using MAT
using Test

@testset "Default section" begin
    # Define the materials:
    materials = [
        Material(29500, 29500, 0.3, 0.3, 11346.154)]

    # Define the nodes:
    nodes = [
        Node(5.00, 1.00, 1),
        Node(5.00, 0.50, 1),
        Node(5.00, 0.00, 1),
        Node(3.75, 0.00, 1),
        Node(2.50, 0.00, 1),
        Node(1.25, 0.00, 1),
        Node(0.00, 0.00, 1),
        Node(0.00, 1.50, 1),
        Node(0.00, 3.00, 1),
        Node(0.00, 4.50, 1),
        Node(0.00, 6.00, 1),
        Node(0.00, 7.50, 1),
        Node(0.00, 9.00, 1),
        Node(1.25, 9.00, 1),
        Node(2.50, 9.00, 1),
        Node(3.75, 9.00, 1),
        Node(5.00, 9.00, 1),
        Node(5.00, 8.50, 1),
        Node(5.00, 8.00, 1)]

    # Define the elements:
    elements = [
        Element(nodes[ 1], nodes[ 2], materials[1], 0.100),
        Element(nodes[ 2], nodes[ 3], materials[1], 0.100),
        Element(nodes[ 3], nodes[ 4], materials[1], 0.100),
        Element(nodes[ 4], nodes[ 5], materials[1], 0.100),
        Element(nodes[ 5], nodes[ 6], materials[1], 0.100),
        Element(nodes[ 6], nodes[ 7], materials[1], 0.100),
        Element(nodes[ 7], nodes[ 8], materials[1], 0.100),
        Element(nodes[ 8], nodes[ 9], materials[1], 0.100),
        Element(nodes[ 9], nodes[10], materials[1], 0.100),
        Element(nodes[10], nodes[11], materials[1], 0.100),
        Element(nodes[11], nodes[12], materials[1], 0.100),
        Element(nodes[12], nodes[13], materials[1], 0.100),
        Element(nodes[13], nodes[14], materials[1], 0.100),
        Element(nodes[14], nodes[15], materials[1], 0.100),
        Element(nodes[15], nodes[16], materials[1], 0.100),
        Element(nodes[16], nodes[17], materials[1], 0.100),
        Element(nodes[17], nodes[18], materials[1], 0.100),
        Element(nodes[18], nodes[19], materials[1], 0.100)]
    
    # Define the section:
    section = Section(materials, nodes, elements)

    # Test the section properties:
    @testset "Section properties" begin
        @test isapprox(   section.A,  2.100000, atol = 1E-6)
        @test isapprox( section.x_c,  1.666667, atol = 1E-6)
        @test isapprox( section.z_c,  4.500000, atol = 1E-6)
        @test isapprox(section.I_xx, 29.541667, atol = 1E-6)
        @test isapprox(section.I_zz,  7.500000, atol = 1E-6)
        @test isapprox(section.I_xz,  0.000000, atol = 1E-6)
    end
end