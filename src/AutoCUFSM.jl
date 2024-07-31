module AutoCUFSM
# --------------------------------------------------
# IMPORT PACKAGES
# --------------------------------------------------
import LinearAlgebra
import StaticArrays
import SparseArrays
import ProgressMeter

# --------------------------------------------------
# DEFINE TYPES
# --------------------------------------------------

# --------------------------------------------------
# EXPORT TYPES AND FUNCTIONS
# --------------------------------------------------
include("Section.jl")
export Material, Node, Element, Section
include("ComputeGeometricProperties.jl")
include("AssembleStiffnessMatrices.jl")
export compute_k_e_local, assemble_K_e_global_sparse
export compute_k_g_local, assemble_K_g_global_sparse
export _transform_k_from_local_to_global
export _compute_undetermined_coefficients
include("ComputeSignatureCurve.jl")
export solve
end
