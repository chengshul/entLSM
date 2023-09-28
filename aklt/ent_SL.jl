using LinearAlgebra

function main1(L::Int, SS::Int)
    d = 2^L
    rho = reshape(collect(Float64,readeach(open("./data/rho/rho_mpo_aklt_S_$(SS/2)_L_$(L).dat", "r"), Float64)), d, d)
    ind = [sum(digits(x,base=2)) for x in 0:2^L-1]
    S = 0
    for m in 0:L
        ind_m = (ind .== m)
        pA = eigvals(Symmetric(rho[ind_m,ind_m]))
        pA = pA[pA .> 1e-12]
        S += -sum(pA.*log.(pA))
    end
    return S
end

function main2(SS::Int)
    L = 16
    S = 0.0
    for m in 0:L
        rho = collect(Float64,readeach(open("./data/rho/rho1_mpo_aklt_S_$(SS/2)_L_$(L)_m_$(m).dat", "r"), Float64))
        d = Int(sqrt(length(rho)))
        rho = reshape(rho, d, d)
        pA = eigvals(Symmetric(rho))
        pA = pA[pA .> 1e-12]
        S += -sum(pA.*log.(pA))
    end
    return S
end

L = parse(Int,ARGS[1])
SS = parse(Int,ARGS[2])
S = 0.0
if L < 16
    S = main1(L,SS)
else
    S = main2(SS)
end

open("./data/SL_mpo1_aklt_S_$(SS/2)_L_$(L).dat", "w") do io
    write(io, [S])
end
