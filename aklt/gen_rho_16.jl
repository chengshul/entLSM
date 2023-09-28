using LinearAlgebra

function main(L::Int, SS::Int)
    d = 2^L
    rho = reshape(collect(Float64,readeach(open("./data/rho/rho0_mpo_aklt_S_$(SS/2)_L_$(L)_i_1.dat", "r"), Float64)), d, d)
    for i in 2:4
        rho += reshape(collect(Float64,readeach(open("./data/rho/rho0_mpo_aklt_S_$(SS/2)_L_$(L)_i_$(i).dat", "r"), Float64)), d, d)
    end
    #open("./data/mpo_M_aklt_S_$(SS/2).dat", "w") do io
    #    write(io, M)
    #end
    open("./data/rho/rho1_mpo_aklt_S_$(SS/2)_L_$(L).dat", "w") do io
        write(io, rho)
    end

    for m in 0:L
        ind = [sum(digits(i, base=2)) for i in 0:2^L-1] .== m
        open("./data/rho/rho1_mpo_aklt_S_$(SS/2)_L_$(L)_m_$(m).dat", "w") do io
            write(io, rho[ind,ind])
        end
    end
end

L = 16 #parse(Int,ARGS[1])
SS = parse(Int,ARGS[2])
main(L,SS)
