function solve(
    section ::Section{NN, NE, NM}, 
    L       ::StaticArrays.SVector{NL, <:Real}, 
    M       ::StaticArrays.SVector{NLT, <:Integer}; 
    NRE     ::Integer = 10) where {NN, NE, NM, NL, NLT}
    elements = section.elements

    Λ = Vector[]
    Φ = Matrix[]
    ProgressMeter.@showprogress desc = "Performing elastic buckling analysis..." for a in L
        K_e = assemble_K_e_global_sparse(Val(NN), elements, a, M)
        K_g = assemble_K_g_global_sparse(Val(NN), elements, a, M)

        λ, ϕ = LinearAlgebra.eigen(Array(K_e), Array(K_g))

        λ = real(λ)
        ϕ = real(ϕ)

        keep_index = findall(x -> x ≥ 0, λ)
        λ = λ[keep_index]
        ϕ = ϕ[:, keep_index]

        sort_index = sortperm(λ)
        λ = λ[sort_index]
        ϕ = ϕ[:, sort_index]

        NRE = min(NRE, length(λ))
        λ = λ[1:NRE]
        ϕ = ϕ[:, 1:NRE]

        for i in 1:NRE
            ϕ[:, i] = ϕ[:, i] / maximum(abs.(ϕ[:, i]))
        end

        push!(Λ, λ)
        push!(Φ, ϕ)
    end

    return Λ, Φ
end

solve(section::Section, L::AbstractVector{<:Real}, M::AbstractVector{<:Integer}) = solve(section, StaticArrays.SVector{length(L)}(L), StaticArrays.SVector{length(M)}(M))