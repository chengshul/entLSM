using LinearAlgebra
using SparseArrays
using Arpack
using Dates

function σ(N::Int)
    d = 2^N
    σz_all = SparseMatrixCSC{Float64, Int}[]
    σp_all = SparseMatrixCSC{Float64, Int}[]
    
    σz0 = sparse([1 0; 0 -1])
    σp0 = sparse([0 1; 0 0])
    Id = sparse([1 0; 0 1])

    for i in 1:N
        σz = spdiagm([1])
        σp = spdiagm([1])
        for j in 1:i-1
            σz = kron(σz, Id)
            σp = kron(σp, Id)
        end
        σz = kron(σz, σz0)
        σp = kron(σp, σp0)
        for j in i+1:N
            σz = kron(σz, Id)
            σp = kron(σp, Id)
        end
        push!(σz_all, σz)
        push!(σp_all, σp)
    end
    return σz_all, σp_all
end


function main(L::Int, SS::Int)
    ind = [sum(digits(i, base=2, pad=L)) for i in 0:2^L-1] .== L÷2
    d = sum(ind)
    ρA = reshape(collect(Float64,readeach(open("./data/rho/rho1_mpo_aklt_S_$(SS/2)_L_$(L)_Sz_0.dat", "r"), Float64)), d, d)

    ES, US = eigen(Symmetric(-ρA))
    ES = -ES
    m = sum(ES .> 1e-14)
    K = -US[:,1:m] * diagm(log.(ES[1:m])) * (US[:,1:m])'
    σz, σp = σ(L)
    psi = US[:,1]
    E0 = psi' * K * psi
    for i in 1:L
        psi = diagm(exp.(1im*π*i*diag(σz[i][ind,ind])/L)) * psi
    end
    E1 = psi' * K * psi
    println([-log(ES[1]),E0,E1])

    open("./data/twist_mpo_aklt_S_$(SS/2)_L_$(L)_Sz_0.dat", "w") do io
        write(io, [E0,E1])
    end

    p1 = psi' * ρA * psi
    println(-log(p1))
    p2 = psi' * ρA * ρA * psi
    println(-log(p2))

    open("./data/twist_mpo_aklt_S_$(SS/2)_L_$(L)_Sz_0_rho.dat", "w") do io
        write(io, [exp(-E0),p1,p2])
    end
end

L = parse(Int,ARGS[1])
SS = parse(Int,ARGS[2])

main(L,SS)
