function compute_local_k_e(a, b, t, E_xx, E_yy, v_xx, v_yy, G_xy, M::StaticArrays.SVector{NLT, <:Integer})::StaticArrays.SMatrix{8 * NLT, 8 * NLT} where {NLT}
    # Preallocate the local elastic stiffness matrix:
    k_e = zeros(8 * NLT, 8 * NLT)

    # Compute the rigidity properties:
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
            k_e_1 = zeros(4, 4) # Membrane elastic stiffness matrix
            k_e_2 = zeros(4, 4) # Flexural elastic stiffness matrix
            k_e_m = zeros(8, 8) # Local elastic stiffness matrix for a single longitudinal term

            # Compute the undetermined coefficients:
            u_i = m_i * π
            u_j = m_j * π
            c_1 = u_i / a
            c_2 = u_j / a
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

            k_e_1 = t * k_e_1

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

            # Define the local elastic stiffness matrix for a single longitudinal term:
            @inbounds k_e_m[1:4, 1:4] = k_e_1
            @inbounds k_e_m[5:8, 5:8] = k_e_2

            # Add:
            @inbounds k_e[(8 * (m_i - 1) + 1):(8 * m_i), (8 * (m_j - 1) + 1):(8 * m_j)] = k_e_m
        end
    end

    # Convert the local elastic stiffness matrix to a static array:
    k_e = StaticArrays.SMatrix{8 * NLT, 8 * NLT}(k_e)

    # Return the result:
    return k_e
end

compute_local_k_e(a, b, t, E_xx, E_yy, v_xx, v_yy, G_xy, M::AbstractVector{<:Integer}) = compute_local_k_e(a, b, t, E_xx, E_yy, v_xx, v_yy, G_xy, StaticArrays.SVector{length(M)}(M))

function compute_local_k_g(a, b, t, σ_i, σ_j, M::StaticArrays.SVector{NLT, <:Integer})::StaticArrays.SMatrix{8 * NLT, 8 * NLT} where {NLT}
    # Preallocate the local geometric stiffness matrix:
    k_g = zeros(8 * NLT, 8 * NLT)

    # Convert stresses to line loads:
    T_i = σ_i * t
    T_j = σ_j * t

    # Assemble the local geometric stiffness matrix:
    for m_i in M
        for m_j in M
            # Preallocate:
            k_g_1 = zeros(4, 4) # Membrane geometric stiffness matrix
            k_g_2 = zeros(4, 4) # Flexural geometric stiffness matrix
            k_g_m = zeros(8, 8) # Local geometric stiffness matrix for a single longitudinal term

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

            # Define the local elastic stiffness matrix for a single longitudinal term:
            @inbounds k_g_m[1:4, 1:4] = k_g_1
            @inbounds k_g_m[5:8, 5:8] = k_g_2

            # Add:
            @inbounds k_g[(8 * (m_i - 1) + 1):(8 * m_i), (8 * (m_j - 1) + 1):(8 * m_j)] = k_g_m
        end
    end

    # Convert the local elastic stiffness matrix to a static array:
    k_g = StaticArrays.SMatrix{8 * NLT, 8 * NLT}(k_g)

    # Return the result:
    return k_g
end

compute_local_k_g(a, b, T_i, T_j, M::AbstractVector{<:Integer}) = compute_local_k_g(a, b, T_i, T_j, StaticArrays.SVector{length(M)}(M))

function assemble_global_K_e(::Val{NN}, elements::StaticArrays.SVector{NE, Element}, a::Real, M::StaticArrays.SVector{NLT, <:Integer})::StaticArrays.SMatrix{4 * NN * NLT, 4 * NN * NLT} where {NN, NE, NLT}
    # Preallocate the global elastic stiffness matrix:
    K_e = zeros(4 * NN * NLT, 4 * NN * NLT)

    # Assemble the global elastic stiffness matrix:
    for element in elements
        # Extract the element information:
        b = element.w
        t = element.t

        # Extract the node information:
        node_i    = element.node_i
        node_j    = element.node_j
        node_i_ID = node_i.ID
        node_j_ID = node_j.ID

        # Extract the material information:
        material = element.material
        E_xx     = material.E_xx
        E_yy     = material.E_yy
        v_xx     = material.v_xx
        v_yy     = material.v_yy
        G_xy     = material.G_xy

        # Compute the local elastic stiffness matrix:
        k_e = compute_local_k_e(a, b, t, E_xx, E_yy, v_xx, v_yy, G_xy, M)

        # Assemble the global elastic stiffness matrix:
        for i in 1:NLT
            for j in 1:NLT
                @inbounds K_e[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_e[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds K_e[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_e[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds K_e[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_e[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds K_e[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_e[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
                @inbounds K_e[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_e[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds K_e[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_e[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds K_e[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_e[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds K_e[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_e[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
                @inbounds K_e[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_e[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds K_e[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_e[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds K_e[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_e[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds K_e[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_e[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
                @inbounds K_e[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_e[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds K_e[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_e[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds K_e[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_e[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds K_e[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_e[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
            end
        end
    end

    # Convert the global elastic stiffness matrix to a static array:
    K_e = StaticArrays.SMatrix{4 * NN * NLT, 4 * NN * NLT}(K_e)

    # Return the result:
    return K_e
end

function assemble_global_K_g(::Val{NN}, elements::StaticArrays.SVector{NE, Element}, a::Real, M::StaticArrays.SVector{NLT, <:Integer})::StaticArrays.SMatrix{4 * NN * NLT, 4 * NN * NLT} where {NN, NE, NLT}
    # Preallocate the global geometric stiffness matrix:
    K_g = zeros(4 * NN * NLT, 4 * NN * NLT)

    # Assemble the global geometric stiffness matrix:
    for element in elements
        # Extract the element information:
        b = element.w
        t = element.t

        # Extract the node information:
        node_i    = element.node_i
        node_j    = element.node_j
        node_i_ID = node_i.ID
        node_j_ID = node_j.ID
        σ_i       = node_i.σ
        σ_j       = node_j.σ

        # Compute the local geometric stiffness matrix:
        k_g = compute_local_k_g(a, b, t, σ_i, σ_j, M)

        # Assemble the global geometric stiffness matrix:
        for i in 1:NLT
            for j in 1:NLT
                @inbounds K_g[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_g[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds K_g[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_g[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds K_g[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_g[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds K_g[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_g[(8 * (i - 1) + 1):(8 * (i - 1) + 2), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
                @inbounds K_g[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_g[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds K_g[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_g[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds K_g[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_g[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds K_g[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_g[(8 * (i - 1) + 3):(8 * (i - 1) + 4), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
                @inbounds K_g[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_g[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds K_g[(4 * NN * (i - 1)          + node_i_ID * 2 - 1):(4 * NN * (i - 1)          + node_i_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_g[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds K_g[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_i_ID * 2)] = k_g[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds K_g[(4 * NN * (i - 1)          + node_j_ID * 2 - 1):(4 * NN * (i - 1)          + node_j_ID * 2), (4 * NN * (j - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (j - 1) + 2 * NN + node_j_ID * 2)] = k_g[(8 * (i - 1) + 5):(8 * (i - 1) + 6), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
                @inbounds K_g[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_g[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 1):(8 * (j - 1) + 2)]
                @inbounds K_g[(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_i_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_g[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 3):(8 * (j - 1) + 4)]
                @inbounds K_g[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1)          + node_i_ID * 2 - 1):(4 * NN * (j - 1)          + node_i_ID * 2)] = k_g[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 5):(8 * (j - 1) + 6)]
                @inbounds K_g[(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2 - 1):(4 * NN * (i - 1) + 2 * NN + node_j_ID * 2), (4 * NN * (j - 1)          + node_j_ID * 2 - 1):(4 * NN * (j - 1)          + node_j_ID * 2)] = k_g[(8 * (i - 1) + 7):(8 * (i - 1) + 8), (8 * (j - 1) + 7):(8 * (j - 1) + 8)]
            end
        end
    end
    
    # Convert the global elastic stiffness matrix to a static array:
    K_g = StaticArrays.SMatrix{4 * NN * NLT, 4 * NN * NLT}(K_g)

    # Return the result:
    return K_g
end

function _compute_undetermined_coefficients(a::T, m_i::Integer, m_j::Integer)::NTuple{5} where {T<:Real}
    F = float(T)

    if m_i == m_j
        I_1 = +a / 2
        I_2 = -(π ^ 2 * m_i ^ 2) / (2 * a)
        I_3 = -(π ^ 2 * m_j ^ 2) / (2 * a)
        I_4 = +(π ^ 4 * m_i ^ 4) / (2 * a ^ 3)
        I_5 = +(π ^ 2 * m_i ^ 2) / (2 * a)
    else
        I_1 = zero(F)
        I_2 = zero(F)
        I_3 = zero(F)
        I_4 = zero(F)
        I_5 = zero(F)
    end

    return I_1, I_2, I_3, I_4, I_5
end