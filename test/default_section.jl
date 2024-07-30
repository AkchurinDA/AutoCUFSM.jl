using AutoCUFSM

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
    Element(nodes,  1,  2, materials, 1, 0.100),
    Element(nodes,  2,  3, materials, 1, 0.100),
    Element(nodes,  3,  4, materials, 1, 0.100),
    Element(nodes,  4,  5, materials, 1, 0.100),
    Element(nodes,  5,  6, materials, 1, 0.100),
    Element(nodes,  6,  7, materials, 1, 0.100),
    Element(nodes,  7,  8, materials, 1, 0.100),
    Element(nodes,  8,  9, materials, 1, 0.100),
    Element(nodes,  9, 10, materials, 1, 0.100),
    Element(nodes, 10, 11, materials, 1, 0.100),
    Element(nodes, 11, 12, materials, 1, 0.100),
    Element(nodes, 12, 13, materials, 1, 0.100),
    Element(nodes, 13, 14, materials, 1, 0.100),
    Element(nodes, 14, 15, materials, 1, 0.100),
    Element(nodes, 15, 16, materials, 1, 0.100),
    Element(nodes, 16, 17, materials, 1, 0.100),
    Element(nodes, 17, 18, materials, 1, 0.100),
    Element(nodes, 18, 19, materials, 1, 0.100)]

# Define the section:
section = Section(materials, nodes, elements)

# Define the lengths:
L = collect(10 .^ range(0, 3, 100))

# Define the number of longitudinal terms:
M = [1]

# Compute the signature curve:
Λ, Φ =  solve(section, L, M)



