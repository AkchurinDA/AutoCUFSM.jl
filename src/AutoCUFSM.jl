module AutoCUFSM
# --------------------------------------------------
# IMPORT PACKAGES
# --------------------------------------------------
import ForwardDiff
import LinearAlgebra
import StaticArrays

# --------------------------------------------------
# DEFINE TYPES
# --------------------------------------------------

# --------------------------------------------------
# EXPORT TYPES AND FUNCTIONS
# --------------------------------------------------
include("Section.jl")
export Section, Node, Element, Material
include("ComputeGeometricProperties.jl")
include("ComputeStiffnesses.jl")
export compute_local_k_e, assemble_global_K_e
export compute_local_k_g, assemble_global_K_g
end
