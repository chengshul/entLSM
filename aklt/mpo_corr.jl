using LinearAlgebra
using WignerSymbols

function rho_f(SS::Int, N::Int)
    A = zeros(2,SS+1,3)
    B = zeros(2,2,3)
    s_1h = [1/2,-1/2]
    s_1 = [1,0,-1]
    s_3h = [3/2,1/2,-1/2,-3/2]
    s_SS = []
    if SS == 1
        s_SS = s_1h
    elseif SS == 3
        s_SS = s_3h
    end

    isy = [0 1; -1 0]
    for i in 1:2
        for j in 1:SS+1
            for k in 1:3
                A[i,j,k] = clebschgordan(1/2,s_1h[i],SS/2,s_SS[j],1,s_1[k])
            end
        end
    end
    for i in 1:2
        for j in 1:2
            for k in 1:3
                B[i,j,k] = clebschgordan(1/2,s_1h[i],1/2,s_1h[j],1,s_1[k])
            end
        end
    end
    B = reshape(B,4,3)
    A = reshape(A,2*(SS+1),3)
    A = reshape(A * B', 4*(SS+1),2)
    A = reshape(A * isy, 2,SS+1,4)
    A = reshape(permutedims(A, (1,3,2)), 8,SS+1)
    A = reshape(A * A', 2,2,2,2,2,2)
    A = reshape(permutedims(A, (4,1,2,5,3,6)), 2,2,4,4)
    A = A * 2 / 3 / ((1+3.0^(1-N))^(1/N))

    return A
end

function corr_f(M::Array{Float64,4}, L::Int)
    M1 = reshape(permutedims(M, (1,3,4,2)), 32,2)
    M1 = reshape(M1 * [1 0; 0 -1], 2,4,4,2)
    M1 = M1[1,:,:,1] + M1[2,:,:,2]
    M2 = M[1,1,:,:] + M[2,2,:,:]
    c = zeros(L-1)
    for i in 1:L-1
        c[i] = tr(M1 * M2^(i-1) * M1 * M2^(L-1-i))
    end
    return c
end

L = parse(Int,ARGS[1])
SS = parse(Int,ARGS[2])
M = rho_f(SS, L)
corr = corr_f(M, L)

open("./data/mpo_corr_aklt_S_$(SS/2)_L_$(L).dat", "w") do io
    write(io, corr)
end
#open("./data/rho/rho_mpo_aklt_S_$(SS/2)_L_$(L).dat", "w") do io
#    write(io, rho)
#end
