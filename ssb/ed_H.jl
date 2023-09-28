using LinearAlgebra
using SparseArrays
using Arpack

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

function ham(L::Int, σz::Vector{SparseMatrixCSC{Float64, Int}},
        σp::Vector{SparseMatrixCSC{Float64, Int}})
    d = 2^L
    H = spzeros(Float64, d, d)
    for i in 1:L
        H += 1/4 * σz[i] * σz[i%L+1]
        H += 1/2 * σp[i] * transpose(σp[i%L+1])
        H += 1/2 * transpose(σp[i]) * σp[i%L+1]
        
        H += 1/8 * σz[i] * σz[(i+1)%L+1]
        H += 1/4 * σp[i] * transpose(σp[(i+1)%L+1])
        H += 1/4 * transpose(σp[i]) * σp[(i+1)%L+1]
    end
    return H
end

function main(L::Int)
    σz, σp = σ(L)
    H = ham(L, σz, σp)
    
    #E, a = eigs(Symmetric(H); nev=40, which=:SR, ritzvec=false)
    E = eigvals(Symmetric(Matrix(H)))
    #println(E)
    open("./data/E_mg_L_$(L).dat", "w") do io
        write(io, E)
    end
end

L = parse(Int,ARGS[1])
main(L)
