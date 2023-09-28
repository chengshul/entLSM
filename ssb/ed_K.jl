using LinearAlgebra
using SparseArrays
using Arpack
using Dates

function Sz_f(N::Int, i::Int)
    L = N÷2
    d = 8^L
    
    Sz_1h_0 = sparse([1/2 0; 0 -1/2])
    Sz_3h_0 = sparse(diagm([3/2,1/2,-1/2,-3/2]))
    Id_1h = sparse(diagm([1.0,1.0]))
    Id_3h = sparse(diagm([1.0,1.0,1.0,1.0]))

    if i <= L
        Sz = spdiagm([1])
        for j in 1:i-1
            Sz = kron(Sz, Id_1h)
        end
        Sz = kron(Sz, Sz_1h_0)
        for j in i+1:L
            Sz = kron(Sz, Id_1h)
        end
        for j in L+1:N
            Sz = kron(Sz, Id_3h)
        end
        return Sz
    else
        Sz = spdiagm([1])
        for j in 1:L
            Sz = kron(Sz, Id_1h)
        end
        for j in 1:i-L-1
            Sz = kron(Sz, Id_3h)
        end
        Sz = kron(Sz, Sz_3h_0)
        for j in i+1-L:L
            Sz = kron(Sz, Id_3h)
        end
        return Sz
    end
end

function Sp_f(N::Int, i::Int)
    L = N÷2
    d = 8^L
    
    Sp_1h_0 = sparse([0 1; 0 0])
    Sp_3h_0 = sparse([0 √3 0 0; 0 0 2 0; 0 0 0 √3; 0 0 0 0]) 
    Id_1h = sparse(diagm([1.0,1.0]))
    Id_3h = sparse(diagm([1.0,1.0,1.0,1.0]))

    if i <= L
        Sp = spdiagm([1])
        for j in 1:i-1
            Sp = kron(Sp, Id_1h)
        end
        Sp = kron(Sp, Sp_1h_0)
        for j in i+1:L
            Sp = kron(Sp, Id_1h)
        end
        for j in L+1:N
            Sp = kron(Sp, Id_3h)
        end
        return Sp
    else
        Sp = spdiagm([1])
        for j in 1:L
            Sp = kron(Sp, Id_1h)
        end
        for j in 1:i-L-1
            Sp = kron(Sp, Id_3h)
        end
        Sp = kron(Sp, Sp_3h_0)
        for j in i+1-L:L
            Sp = kron(Sp, Id_3h)
        end
        return Sp
    end
end

function ham(N::Int, J1::Float64, J2::Float64, D::Float64)
    L = N÷2
    d = 8^L
    H = spzeros(Float64, d, d)
    for i in 1:L
        H += J1 * (Sz_f(N,i) * Sz_f(N,i%L+1)
                   +Sp_f(N,i) * transpose(Sp_f(N,i%L+1))/2
                   +transpose(Sp_f(N,i)) * Sp_f(N,i%L+1)/2)
        H += J1/2 * (Sz_f(N,i) * Sz_f(N,(i+1)%L+1)
                     +Sp_f(N,i) * transpose(Sp_f(N,(i+1)%L+1))/2
                     +transpose(Sp_f(N,i)) * Sp_f(N,(i+1)%L+1)/2)
        H += J2 * (Sz_f(N,i) * Sz_f(N,i+L)
                   +Sp_f(N,i) * transpose(Sp_f(N,i+L))/2
                   +transpose(Sp_f(N,i)) * Sp_f(N,i+L)/2)
        H += D * (Sz_f(N,i)+Sz_f(N,i+L))^2
    end
    return H
end

function rhoA(N::Int, ψ::Vector{Float64}, LA::Int)
    LB = N÷2*3 - LA
    
    ρA = zeros(Float64, 2^LA, 2^LA)
    for i in 0:2^(N÷2*3)-1
        if i % 100000 == 0
            println(i)
            println(Dates.format(now(), "yyyy u d HH:MM:SS"))
            flush(stdout)
        end
        i_b = digits(i, base=2, pad=N÷2*3)
        i_0 = i % (2^LB)
        for j in 0:2^LA-1
            j_b = digits(j, base=2, pad=LA)
            ρA[evalpoly(2, i_b[LB+1:N÷2*3])+1, evalpoly(2, j_b)+1] += ψ[i+1,1] * ψ[j*2^LB+i_0+1,1]
        end
    end
    println(issymmetric(ρA))
    return ρA
end

function main(N::Int, J1::Float64, J2::Float64, D::Float64)
    H = ham(N, J1, J2, D)
    E, ψ = eigs(Symmetric(H); nev=10, which=:SR)
    println(E)
    
    #=z = 0.0
    for i in 1:N
        z = z .+ Sz_f(N,i) * ψ
    end
    z = ψ' * z
    println(z)=#

    open("./data/E_ssb_3h_z_N_$(N)_J1_$(J1)_J2_$(J2)_D_$(D).dat", "w") do io
        write(io, E)
    end
    
    l = N÷2
    ρA = rhoA(N, ψ[:,1], l)
    open("./data/rho_ssb_3h_z_N_$(N)_J1_$(J1)_J2_$(J2)_D_$(D)_l_$(l).dat", "w") do io
        write(io, ρA)
    end
    ES, US = eigen(Symmetric(ρA))
    println(ES[end:-1:end-10])
    open("./data/ES_ssb_3h_z_N_$(N)_J1_$(J1)_J2_$(J2)_D_$(D)_l_$(l).dat", "w") do io
        write(io, ES)
    end
    println(-log.(ES)[end:-1:end-10])
end

L = parse(Int,ARGS[1])
N = 2*L
J1 = parse(Float64,ARGS[2])
J2 = parse(Float64,ARGS[3])
D = parse(Float64,ARGS[4])

main(N, J1, J2, D)
