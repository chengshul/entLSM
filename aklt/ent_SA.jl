using LinearAlgebra
using WignerSymbols

function rho_f(N::Int, SS::Int)
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

    B = zeros(4,4)
    for i in 1:2
        B += A[i,i,:,:]
    end

    #=rho = zeros(2^N,2^N)
    for i in 1:4
        rho0 = A[:,:,i,:]
        for j in 1:N-2
            rho0 = reshape(rho0,4^j,4)*reshape(permutedims(A,[3,1,2,4]),4,16)
            rho0 = permutedims(reshape(rho0,2^j,2^j,2,2,4),[1,3,2,4,5])
            rho0 = reshape(rho0,2^(j+1),2^(j+1),4)
        end
        rho0 = reshape(rho0,4^(N-1),4)*reshape(permutedims(A[:,:,:,i],[3,1,2]),4,4)
        rho0 = permutedims(reshape(rho0,2^(N-1),2^(N-1),2,2),[1,3,2,4])
        rho += reshape(rho0,2^N,2^N)
    end=#
    
    SA_l = Float64[]
    for l in 1:N-1
        rhoA = zeros(2^l,2^l)
        for i in 1:4
            rho0 = A[:,:,i,:]
            for j in 1:l-1
                rho0 = reshape(rho0,4^j,4)*reshape(permutedims(A,[3,1,2,4]),4,16)
                rho0 = permutedims(reshape(rho0,2^j,2^j,2,2,4),[1,3,2,4,5])
                rho0 = reshape(rho0,2^(j+1),2^(j+1),4)
            end
            rho0 = reshape(rho0,4^l,4)
            for j in 1:N-l-1
                rho0 = rho0 * B
            end
            rhoA += reshape(rho0 * B[:,i],2^l,2^l)
        end
        #println(tr(rhoA))
        
        ind = [sum(digits(x,base=2)) for x in 0:2^l-1]
        SA = 0
        for m in 0:l
            ind_m = (ind .== m)
            pA = eigvals(Symmetric(rhoA[ind_m,ind_m]))
            pA = pA[pA .> 1e-12]
            SA += -sum(pA.*log.(pA))
        end
        push!(SA_l,SA)
    end

    SB_l = Float64[]
    for l in 1:(N-1)÷3
        rhoB = zeros(2^(2*l),2^(2*l))
        for i in 1:4
            rho0 = A[:,:,i,:]
            for j in 1:l-1
                rho0 = reshape(rho0,4^j,4)*reshape(permutedims(A,[3,1,2,4]),4,16)
                rho0 = permutedims(reshape(rho0,2^j,2^j,2,2,4),[1,3,2,4,5])
                rho0 = reshape(rho0,2^(j+1),2^(j+1),4)
            end
            rho0 = reshape(rho0,4^l,4)
            for j in 1:l
                rho0 = rho0 * B
            end
            for j in l:2*l-1
                rho0 = reshape(rho0,4^j,4)*reshape(permutedims(A,[3,1,2,4]),4,16)
                rho0 = permutedims(reshape(rho0,2^j,2^j,2,2,4),[1,3,2,4,5])
                rho0 = reshape(rho0,2^(j+1),2^(j+1),4)
            end
            rho0 = reshape(rho0,4^(2*l),4)
            for j in 1:N-3*l-1
                rho0 = rho0 * B
            end
            rhoB += reshape(rho0 * B[:,i], 2^(2*l),2^(2*l))
        end
        #println(rhoB == rhoB')
        #println(tr(rhoB))
        ind = [sum(digits(x,base=2)) for x in 0:2^(2*l)-1]
        SB = 0
        for m in 0:2*l
            ind_m = (ind .== m)
            pB = eigvals(Symmetric(rhoB[ind_m,ind_m]))
            pB = pB[pB .> 1e-12]
            SB += -sum(pB.*log.(pB))
        end
        push!(SB_l,SB)
    end
    return A,SA_l,SB_l
end

L = parse(Int,ARGS[1])
SS = parse(Int,ARGS[2])
M,SA,SB = rho_f(L,SS)

open("./data/SA_mpo1_aklt_S_$(SS/2)_L_$(L).dat", "w") do io
    write(io, SA)
end
open("./data/S2A_mpo1_aklt_S_$(SS/2)_L_$(L).dat", "w") do io
    write(io, SB)
end
#open("./data/mpo_M_aklt_S_$(SS/2).dat", "w") do io
#    write(io, M)
#end
#open("./data/rho/rho_mpo_aklt_S_$(SS/2)_L_$(L).dat", "w") do io
#    write(io, rho)
#end
