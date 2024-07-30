function compute_k_e_local(
    a   ::GeometricPropertyType, 
    b   ::GeometricPropertyType, 
    t   ::GeometricPropertyType, 
    E_xx::MaterialPropertyType,
    E_yy::MaterialPropertyType,
    v_xx::MaterialPropertyType,
    v_yy::MaterialPropertyType,
    G_xy::MaterialPropertyType,
    M   ::StaticArrays.SVector{NLT, <:Integer})::StaticArrays.SMatrix{8 * NLT, 8 * NLT, GeometricPropertyType} where {
        GeometricPropertyType   <: Real, 
        MaterialPropertyType    <: Real, 
        NLT}
    # Preallocate the local elastic stiffness matrix:
    k_e = zeros(GeometricPropertyType, 8 * NLT, 8 * NLT)

    # Compute element's rigidities:
    E_11 = E_xx / (1 - v_xx * v_yy)
    E_22 = E_yy / (1 - v_xx * v_yy)
    D_xx = (E_xx * t ^ 3) / (12 * (1 - v_xx * v_yy))
    D_yy = (E_yy * t ^ 3) / (12 * (1 - v_xx * v_yy))
    D_11 = (E_yy * v_xx * t ^ 3) / (12 * (1 - v_xx * v_yy))
    D_xy = (G_xy * t ^ 3) / 12

    # Compute the local elastic stiffness matrix:
    for m_i in M
        for m_j in M
            # Preallocate:
            k_e_1 = zeros(GeometricPropertyType, 4, 4) # Membrane elastic stiffness matrix
            k_e_2 = zeros(GeometricPropertyType, 4, 4) # Flexural elastic stiffness matrix

            # Compute the undetermined coefficients:
            c_1 = m_i * π / a
            c_2 = m_j * π / a
            I_1, I_2, I_3, I_4, I_5 = _compute_undetermined_coefficients(a, m_i, m_j)

            # Define the membrane elastic stiffness matrix:
            @inbounds k_e_1[1, 1] = +(E_11            * I_1) / (1 * b            ) + (G_xy * b * I_5) / (3                )
            @inbounds k_e_1[1, 2] = -(E_22 * v_xx     * I_3) / (2           * c_2) - (G_xy *     I_5) / (2           * c_2)
            @inbounds k_e_1[1, 3] = -(E_11            * I_1) / (1 * b            ) + (G_xy * b * I_5) / (6                )
            @inbounds k_e_1[1, 4] = -(E_22 * v_xx     * I_3) / (2           * c_2) + (G_xy *     I_5) / (2           * c_2)
            @inbounds k_e_1[2, 1] = -(E_22 * v_xx     * I_2) / (2     * c_1      ) - (G_xy *     I_5) / (2     * c_1      )
            @inbounds k_e_1[2, 2] = +(E_22        * b * I_4) / (3     * c_1 * c_2) + (G_xy *     I_5) / (1 * b * c_1 * c_2)
            @inbounds k_e_1[2, 3] = +(E_22 * v_xx     * I_2) / (2     * c_1      ) - (G_xy *     I_5) / (2     * c_1      )
            @inbounds k_e_1[2, 4] = +(E_22        * b * I_4) / (6     * c_1 * c_2) - (G_xy *     I_5) / (1 * b * c_1 * c_2)
            @inbounds k_e_1[3, 1] = -(E_11            * I_1) / (b                ) + (G_xy * b * I_5) / (6                )
            @inbounds k_e_1[3, 2] = +(E_22 * v_xx     * I_3) / (2           * c_2) - (G_xy *     I_5) / (2           * c_2)
            @inbounds k_e_1[3, 3] = k_e_1[1, 1]
            @inbounds k_e_1[3, 4] = +(E_22 * v_xx     * I_3) / (2           * c_2) + (G_xy *     I_5) / (2           * c_2)
            @inbounds k_e_1[4, 1] = -(E_22 * v_xx     * I_2) / (2     * c_1      ) + (G_xy *     I_5) / (2     * c_1      )
            @inbounds k_e_1[4, 2] = k_e_1[2, 4]
            @inbounds k_e_1[4, 3] = +(E_22 * v_xx     * I_2) / (2     * c_1      ) + (G_xy *     I_5) / (2     * c_1      )
            @inbounds k_e_1[4, 4] = k_e_1[2, 2]
            @inbounds k_e_1 = t * k_e_1

            # Define the flexural elastic stiffness matrix:
            @inbounds k_e_2[1, 1] = (+5040         * D_xx * I_1 - 504 * b ^ 2 * D_11 * I_2 - 504 * b ^ 2 * D_11 * I_3 + 156 * b ^ 4 * D_yy * I_4 + 2016 * b ^ 2 * D_xy * I_5) / (420 * b ^ 3)
            @inbounds k_e_2[1, 2] = (+2520 * b     * D_xx * I_1 - 462 * b ^ 3 * D_11 * I_2 -  42 * b ^ 3 * D_11 * I_3 +  22 * b ^ 5 * D_yy * I_4 +  168 * b ^ 3 * D_xy * I_5) / (420 * b ^ 3)
            @inbounds k_e_2[1, 3] = (-5040         * D_xx * I_1 + 504 * b ^ 2 * D_11 * I_2 + 504 * b ^ 2 * D_11 * I_3 +  54 * b ^ 4 * D_yy * I_4 - 2016 * b ^ 2 * D_xy * I_5) / (420 * b ^ 3)
            @inbounds k_e_2[1, 4] = (+2520 * b     * D_xx * I_1 -  42 * b ^ 3 * D_11 * I_2 -  42 * b ^ 3 * D_11 * I_3 -  13 * b ^ 5 * D_yy * I_4 +  168 * b ^ 3 * D_xy * I_5) / (420 * b ^ 3)
            @inbounds k_e_2[2, 1] = (+2520 * b     * D_xx * I_1 - 462 * b ^ 3 * D_11 * I_3 -  42 * b ^ 3 * D_11 * I_2 +  22 * b ^ 5 * D_yy * I_4 +  168 * b ^ 3 * D_xy * I_5) / (420 * b ^ 3)
            @inbounds k_e_2[2, 2] = (+1680 * b ^ 2 * D_xx * I_1 -  56 * b ^ 4 * D_11 * I_2 -  56 * b ^ 4 * D_11 * I_3 +   4 * b ^ 6 * D_yy * I_4 +  224 * b ^ 4 * D_xy * I_5) / (420 * b ^ 3)
            @inbounds k_e_2[2, 3] = (-2520 * b     * D_xx * I_1 +  42 * b ^ 3 * D_11 * I_2 +  42 * b ^ 3 * D_11 * I_3 +  13 * b ^ 5 * D_yy * I_4 -  168 * b ^ 3 * D_xy * I_5) / (420 * b ^ 3)
            @inbounds k_e_2[2, 4] = ( +840 * b ^ 2 * D_xx * I_1 +  14 * b ^ 4 * D_11 * I_2 +  14 * b ^ 4 * D_11 * I_3 -   3 * b ^ 6 * D_yy * I_4 -   56 * b ^ 4 * D_xy * I_5) / (420 * b ^ 3)
            @inbounds k_e_2[3, 1] = k_e_2[1, 3]
            @inbounds k_e_2[3, 2] = k_e_2[2, 3]
            @inbounds k_e_2[3, 3] = k_e_2[1, 1]
            @inbounds k_e_2[3, 4] = (-2520 * b     * D_xx * I_1 + 462 * b ^ 3 * D_11 * I_2 +  42 * b ^ 3 * D_11 * I_3 -  22 * b ^ 5 * D_yy * I_4 -  168 * b ^ 3 * D_xy * I_5) / (420 * b ^ 3)
            @inbounds k_e_2[4, 1] = k_e_2[1, 4]
            @inbounds k_e_2[4, 2] = k_e_2[2, 4]
            @inbounds k_e_2[4, 3] = (-2520 * b     * D_xx * I_1 + 462 * b ^ 3 * D_11 * I_3 +  42 * b ^ 3 * D_11 * I_2 -  22 * b ^ 5 * D_yy * I_4 -  168 * b ^ 3 * D_xy * I_5) / (420 * b ^ 3)
            @inbounds k_e_2[4, 4] = k_e_2[2, 2]

            # Add:
            @inbounds k_e[(8 * (m_i - 1) + 1):(8 * (m_i - 1) + 4), (8 * (m_j - 1) + 1):(8 * (m_j - 1) + 4)] = k_e_1
            @inbounds k_e[(8 * (m_i - 1) + 5):(8 * (m_i - 1) + 8), (8 * (m_j - 1) + 5):(8 * (m_j - 1) + 8)] = k_e_2
        end
    end

    # Convert to a static matrix:
    k_e = StaticArrays.SMatrix{8 * NLT, 8 * NLT, GeometricPropertyType}(k_e)

    # Return the result:
    return k_e
end

function compute_k_g_local(
    a   ::GeometricPropertyType, 
    b   ::GeometricPropertyType, 
    t   ::GeometricPropertyType, 
    σ_i ::StressType, 
    σ_j ::StressType, 
    M   ::StaticArrays.SVector{NLT, <:Integer})::StaticArrays.SMatrix{8 * NLT, 8 * NLT, GeometricPropertyType} where {
        GeometricPropertyType   <: Real,
        StressType              <: Real,
        NLT}
    # Preallocate the local geometric stiffness matrix:
    k_g = zeros(GeometricPropertyType, 8 * NLT, 8 * NLT)

    # Convert stresses to line loads:
    T_i = σ_i * t
    T_j = σ_j * t

    # Assemble the local geometric stiffness matrix:
    for m_i in M
        for m_j in M
            # Preallocate:
            k_g_1 = zeros(GeometricPropertyType, 4, 4) # Membrane geometric stiffness matrix
            k_g_2 = zeros(GeometricPropertyType, 4, 4) # Flexural geometric stiffness matrix

            # Compute the undetermined coefficients:
            u_i = m_i * π
            u_j = m_j * π
            _, _, _, I_4, I_5 = _compute_undetermined_coefficients(a, m_i, m_j)

            # Define the membrane geometric stiffness matrix:
            @inbounds k_g_1[1, 1] = +b         * (3 * T_i + 1 * T_j) * I_5 / (12            )
            @inbounds k_g_1[1, 3] = +b         * (1 * T_i + 1 * T_j) * I_5 / (12            )
            @inbounds k_g_1[2, 2] = +b * a ^ 2 * (3 * T_i + 1 * T_j) * I_4 / (12 * u_i * u_j)
            @inbounds k_g_1[2, 4] = +b * a ^ 2 * (1 * T_i + 1 * T_j) * I_4 / (12 * u_i * u_j)
            @inbounds k_g_1[3, 1] = k_g_1[1, 3]
            @inbounds k_g_1[3, 3] = +b         * (1 * T_i + 3 * T_j) * I_5 / (12            )
            @inbounds k_g_1[4, 2] = k_g_1[2, 4]
            @inbounds k_g_1[4, 4] = +b * a ^ 2 * (1 * T_i + 3 * T_j) * I_4 / (12 * u_i * u_j)

            # Define the flexural geometric stiffness matrix:
            @inbounds k_g_2[1, 1] = +(10 * T_i +  3 * T_j) * b     * I_5 /  35
            @inbounds k_g_2[1, 2] = +(15 * T_i +  7 * T_j) * b ^ 2 * I_5 / 420
            @inbounds k_g_2[1, 3] = +( 9 * T_i +  9 * T_j) * b     * I_5 / 140
            @inbounds k_g_2[1, 4] = -( 7 * T_i +  6 * T_j) * b ^ 2 * I_5 / 420
            @inbounds k_g_2[2, 1] = k_g_2[1, 2]
            @inbounds k_g_2[2, 2] = +( 5 * T_i +  3 * T_j) * b ^ 3 * I_5 / 840
            @inbounds k_g_2[2, 3] = +( 6 * T_i +  7 * T_j) * b ^ 2 * I_5 / 420
            @inbounds k_g_2[2, 4] = -( 1 * T_i +  1 * T_j) * b ^ 3 * I_5 / 280
            @inbounds k_g_2[3, 1] = k_g_2[1, 3]
            @inbounds k_g_2[3, 2] = k_g_2[2, 3]
            @inbounds k_g_2[3, 3] = +( 3 * T_i + 10 * T_j) * b     * I_5 /  35
            @inbounds k_g_2[3, 4] = -( 7 * T_i + 15 * T_j) * b ^ 2 * I_5 / 420
            @inbounds k_g_2[4, 1] = k_g_2[1, 4]
            @inbounds k_g_2[4, 2] = k_g_2[2, 4]
            @inbounds k_g_2[4, 3] = k_g_2[3, 4]
            @inbounds k_g_2[4, 4] = +( 3 * T_i +  5 * T_j) * b ^ 3 * I_5 / 840

            # Add:
            @inbounds k_g[(8 * (m_i - 1) + 1):(8 * (m_i - 1) + 4), (8 * (m_j - 1) + 1):(8 * (m_j - 1) + 4)] = k_g_1
            @inbounds k_g[(8 * (m_i - 1) + 5):(8 * (m_i - 1) + 8), (8 * (m_j - 1) + 5):(8 * (m_j - 1) + 8)] = k_g_2
        end
    end

    # Convert to a static matrix:
    k_g = StaticArrays.SMatrix{8 * NLT, 8 * NLT, GeometricPropertyType}(k_g)

    # Return the result:
    return k_g
end

function assemble_K_e_global_sparse(::Val{NN}, 
    elements::StaticArrays.SVector{NE, <:Element},
    a::Real, 
    M::StaticArrays.SVector{NLT, <:Integer})::SparseArrays.SparseMatrixCSC where {NN, NE, NLT}
    # Preallocate the global elastic stiffness matrix:
    K_e = SparseArrays.spzeros(4 * NN * NLT, 4 * NN * NLT)

    # Assemble the global elastic stiffness matrix:
    for (i, element) in enumerate(elements)
        # Extract the element information:
        b = element.w
        θ = element.θ
        t = element.t

        # Extract the node information:
        node_i_ID = element.node_i_ID
        node_j_ID = element.node_j_ID

        # Extract the material information:
        material = element.material
        E_xx     = material.E_xx
        E_yy     = material.E_yy
        v_xx     = material.v_xx
        v_yy     = material.v_yy
        G_xy     = material.G_xy

        # Compute the local elastic stiffness matrix:
        geometric_properties = promote(a, b, t, θ)
        k_e_local  = compute_k_e_local(geometric_properties[1:3]..., E_xx, E_yy, v_xx, v_yy, G_xy, M)
        k_e_global = _transform_k_e_from_local_to_global(Val(NLT), k_e_local, geometric_properties[4])

        # Assemble the global elastic stiffness matrix:
        K_e_temp = SparseArrays.spzeros(4 * NN * NLT, 4 * NN * NLT)
        for i in 1:NLT
            for j in 1:NLT
                @inbounds k_11 = k_e_global[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds k_12 = k_e_global[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds k_13 = k_e_global[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds k_14 = k_e_global[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
                @inbounds k_21 = k_e_global[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds k_22 = k_e_global[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds k_23 = k_e_global[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds k_24 = k_e_global[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
                @inbounds k_31 = k_e_global[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds k_32 = k_e_global[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds k_33 = k_e_global[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds k_34 = k_e_global[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
                @inbounds k_41 = k_e_global[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds k_42 = k_e_global[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds k_43 = k_e_global[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds k_44 = k_e_global[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]

                @inbounds K_e_temp[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_11
                @inbounds K_e_temp[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_12
                @inbounds K_e_temp[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_21
                @inbounds K_e_temp[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_22
                @inbounds K_e_temp[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_33
                @inbounds K_e_temp[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_34
                @inbounds K_e_temp[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_43
                @inbounds K_e_temp[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_44
                @inbounds K_e_temp[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_13
                @inbounds K_e_temp[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_14
                @inbounds K_e_temp[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_23
                @inbounds K_e_temp[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_24
                @inbounds K_e_temp[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_31
                @inbounds K_e_temp[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_32
                @inbounds K_e_temp[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_41
                @inbounds K_e_temp[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_42
            end
        end

        K_e += K_e_temp
    end

    # Return the result:
    return K_e
end

function assemble_K_g_global_sparse(::Val{NN}, 
    elements::StaticArrays.SVector{NE, <:Element},
    a::Real, 
    M::StaticArrays.SVector{NLT, <:Integer})::SparseArrays.SparseMatrixCSC where {NN, NE, NLT}
    # Preallocate the global geometric stiffness matrix:
    K_g = SparseArrays.spzeros(4 * NN * NLT, 4 * NN * NLT)

    # Assemble the global geometric stiffness matrix:
    for (i, element) in enumerate(elements)
        # Extract the element information:
        b = element.w
        θ = element.θ
        t = element.t

        # Extract the node information:
        node_i_ID = element.node_i_ID
        node_j_ID = element.node_j_ID
        σ_i       = element.node_i.σ
        σ_j       = element.node_j.σ

        # Compute the element geometric stiffness matrix:
        geometric_properties = promote(a, b, t, θ)
        k_g_local  = compute_k_g_local(geometric_properties[1:3]..., σ_i, σ_j, M)
        k_g_global = _transform_k_g_from_local_to_global(Val(NLT), k_g_local, geometric_properties[4])

        # Assemble the global geometric stiffness matrix:
        K_g_temp = SparseArrays.spzeros(4 * NN * NLT, 4 * NN * NLT)
        for i in 1:NLT
            for j in 1:NLT
                @inbounds k_11 = k_g_global[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds k_12 = k_g_global[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds k_13 = k_g_global[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds k_14 = k_g_global[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
                @inbounds k_21 = k_g_global[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds k_22 = k_g_global[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds k_23 = k_g_global[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds k_24 = k_g_global[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
                @inbounds k_31 = k_g_global[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds k_32 = k_g_global[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds k_33 = k_g_global[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds k_34 = k_g_global[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
                @inbounds k_41 = k_g_global[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds k_42 = k_g_global[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds k_43 = k_g_global[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds k_44 = k_g_global[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]

                @inbounds K_g_temp[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_11
                @inbounds K_g_temp[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_12
                @inbounds K_g_temp[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_21
                @inbounds K_g_temp[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_22
                @inbounds K_g_temp[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_33
                @inbounds K_g_temp[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_34
                @inbounds K_g_temp[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_43
                @inbounds K_g_temp[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_44
                @inbounds K_g_temp[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_13
                @inbounds K_g_temp[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_14
                @inbounds K_g_temp[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_23
                @inbounds K_g_temp[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_24
                @inbounds K_g_temp[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_31
                @inbounds K_g_temp[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_32
                @inbounds K_g_temp[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_41
                @inbounds K_g_temp[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_42
            end
        end
        
        K_g += K_g_temp
    end

    # Return the result:
    return K_g
end

function _compute_undetermined_coefficients(a::T, m_i::Integer, m_j::Integer)::NTuple{5} where {T<:Real}
    if m_i == m_j
        I_1 = +a / 2
        I_2 = -(π ^ 2 * m_i ^ 2) / (2 * a)
        I_3 = -(π ^ 2 * m_j ^ 2) / (2 * a)
        I_4 = +(π ^ 4 * m_i ^ 4) / (2 * a ^ 3)
        I_5 = +(π ^ 2 * m_i ^ 2) / (2 * a)
    else
        I_1 = zero(T)
        I_2 = zero(T)
        I_3 = zero(T)
        I_4 = zero(T)
        I_5 = zero(T)
    end

    return I_1, I_2, I_3, I_4, I_5
end

function _transform_k_e_from_local_to_global(::Val{NLT}, k_e_local::StaticArrays.SMatrix, θ::GeometricPropertyType)::StaticArrays.SMatrix{8 * NLT, 8 * NLT} where {NLT, GeometricPropertyType <: Real}
    r = SparseArrays.spzeros(GeometricPropertyType, 8, 8)
    r[1, 1] = +cos(θ)
    r[1, 5] = -sin(θ)
    r[2, 2] = 1
    r[3, 3] = +cos(θ)
    r[3, 7] = -sin(θ)
    r[4, 4] = 1
    r[5, 1] = +sin(θ)
    r[5, 5] = +cos(θ)
    r[6, 6] = 1
    r[7, 3] = +sin(θ)
    r[7, 7] = +cos(θ)
    r[8, 8] = 1
    
    R = SparseArrays.spzeros(GeometricPropertyType, 8 * NLT, 8 * NLT)
    for i in 1:NLT
        R[8 * (i - 1) + 1:8 * i, 8 * (i - 1) + 1:8 * i] = r
    end

    k_e_global = LinearAlgebra.transpose(R) * k_e_local * R

    k_e_global = StaticArrays.SMatrix{8 * NLT, 8 * NLT, GeometricPropertyType}(k_e_global)

    return k_e_global
end

function _transform_k_g_from_local_to_global(::Val{NLT}, k_g_local::StaticArrays.SMatrix, θ::GeometricPropertyType)::StaticArrays.SMatrix{8 * NLT, 8 * NLT} where {NLT, GeometricPropertyType <: Real}
    r = SparseArrays.spzeros(GeometricPropertyType, 8, 8)
    r[1, 1] = +cos(θ)
    r[1, 5] = -sin(θ)
    r[2, 2] = 1
    r[3, 3] = +cos(θ)
    r[3, 7] = -sin(θ)
    r[4, 4] = 1
    r[5, 1] = +sin(θ)
    r[5, 5] = +cos(θ)
    r[6, 6] = 1
    r[7, 3] = +sin(θ)
    r[7, 7] = +cos(θ)
    r[8, 8] = 1
    
    R = SparseArrays.spzeros(GeometricPropertyType, 8 * NLT, 8 * NLT)
    for i in 1:NLT
        R[8 * (i - 1) + 1:8 * i, 8 * (i - 1) + 1:8 * i] = r
    end

    k_g_global = R * k_g_local * R

    k_g_global = StaticArrays.SMatrix{8 * NLT, 8 * NLT, GeometricPropertyType}(k_g_global)

    return k_g_global
end