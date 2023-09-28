using LinearAlgebra
using SparseArrays
using Arpack
using Dates

function σ(L::Int)
    d = 2^L
    σz_all = SparseMatrixCSC{Float64, Int}[]
    σp_all = SparseMatrixCSC{Float64, Int}[]
    
    σz0 = sparse([1 0; 0 -1])
    σp0 = sparse([0 1; 0 0])
    Id = sparse([1 0; 0 1])

    for i in 1:L
        σz = spdiagm([1])
        σp = spdiagm([1])
        for j in 1:i-1
            σz = kron(σz, Id)
            σp = kron(σp, Id)
        end
        σz = kron(σz, σz0)
        σp = kron(σp, σp0)
        for j in i+1:L
            σz = kron(σz, Id)
            σp = kron(σp, Id)
        end
        push!(σz_all, σz)
        push!(σp_all, σp)
    end
    return σz_all, σp_all
end

function S2_gen(L::Int, 
        σz::Vector{SparseMatrixCSC{Float64, Int}},
        σp::Vector{SparseMatrixCSC{Float64, Int}})
    d = 2^L
    S2 = spzeros(Float64, d, d)
    Sp = spzeros(Float64, d, d)
    Sz = spzeros(Float64, d, d)
    for i in 1:L
        Sp += σp[i]
        Sz += σz[i]
    end
    S2 = 1/4*Sz*Sz + 1/2*Sp*transpose(Sp) + 1/2*transpose(Sp)*Sp
    return S2
end

function T_gen(L::Int)
    d = 2^L
    T = zeros(Float64, d, d)
    for i in 0:d-1
        j = evalpoly(2, circshift(digits(i, base=2, pad=L), 1))
        T[i+1,j+1] = 1
    end
    return T
end

function main(L::Int, SS::Int)
    d = 2^L
    σz, σp = σ(L)
    S2 = S2_gen(L, σz, σp)
    T = T_gen(L)
    rho = reshape(collect(Float64,readeach(open("./data/rho/rho_mpo_aklt_S_$(SS/2)_L_$(L).dat", "r"), Float64)), d, d)
    E, ψ = eigen(Symmetric(-rho))
    E = -log.(abs.(E))
    
    E_c = [E[1]]
    i_l = [1]
    c = 1
    i = 2
    while c <= 40
        if E[i] - E_c[c] > 1e-8
            append!(E_c,E[i])
            append!(i_l,i)
            c += 1
        end
        i += 1
    end

    #println(E[1:100])
    #println(E_c)
    #println(i_l)

    ES_1 = Float64[]
    ES_k = Int[]
    ES_Stot = Int[]

    for i in 1:40
        U = ψ[:,i_l[i] : i_l[i+1]-1]
        TU = transpose(U) * T * U
        tu = mod.(round.(Int,real(log.(Complex.(eigvals(TU)))/2im/π*L)),L)
        S2U = transpose(U) * S2 * U
        s2u = eigvals(Symmetric(S2U))
        s2u = round.(Int,(sqrt.(1 .+ 4*s2u) .- 1) / 2)
        #println("k",tu)
        #println("Stot",s2u)

        if 2*s2u[1]+1 == length(tu)
            append!(ES_1, E_c[i])
            append!(ES_k, tu[1])
            append!(ES_Stot, s2u[1])
        elseif 2*(2*s2u[1]+1) == length(tu)
            append!(ES_1, E_c[i])
            append!(ES_1, E_c[i])
            append!(ES_k, minimum(tu))
            append!(ES_k, maximum(tu))
            append!(ES_Stot, s2u[1])
            append!(ES_Stot, s2u[1])
        else
            println("higher degeneracy!")
        end
    end

    open("./data/ES_1_mpo_aklt_S_$(SS/2)_L_$(L).dat", "w") do io
        write(io, ES_1)
    end
    open("./data/ES_Stot_mpo_aklt_S_$(SS/2)_L_$(L).dat", "w") do io
        write(io, ES_Stot)
    end
    open("./data/ES_k_mpo_aklt_S_$(SS/2)_L_$(L).dat", "w") do io
        write(io, ES_k)
    end
end

L = parse(Int,ARGS[1])
SS = parse(Int,ARGS[2])

main(L, SS)
