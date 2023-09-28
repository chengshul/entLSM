using ITensors

ITensors.op(::OpName"sz") = [1 0;0 -1]

function main(N::Int, SS::Int)
    ρ = MPO(N)
    coe = 1.0
    if SS == 1
        coe = 0.8
    else
        coe = 0.6
    end
    T0 = reshape(collect(Float64,readeach(open("./data/mpo_M_aklt_S_$(SS/2)_L_$(N).dat", "r"), Float64)), 2,2,4,4) / coe
    T0 = permutedims(T0, [3,4,1,2])

    T00 = zeros(16,2,2)
    T01 = zeros(16,16,2,2)
    T02 = zeros(16,2,2)
    for i in 0:3
        T00[i*4+1:(i+1)*4,:,:] = T0[i+1,:,:,:]
        T01[i*4+1:(i+1)*4,i*4+1:(i+1)*4,:,:] = T0
        T02[i*4+1:(i+1)*4,:,:] = T0[:,i+1,:,:]
    end
   
    ind_p = Index[]
    ind_b = Index[]
    for i in 1:N
        push!(ind_p, Index(QN("m",-1)=>1,QN("m",1)=>1;tags="p$(i)"))
    end
    for i in 1:N-1
        push!(ind_b, Index(QN("m",0)=>1,
                           QN("m",2)=>1,
                           QN("m",-2)=>1,
                           QN("m",0)=>1,
                           QN("m",-2)=>1,
                           QN("m",0)=>1,
                           QN("m",-4)=>1,
                           QN("m",-2)=>1,
                           QN("m",2)=>1,
                           QN("m",4)=>1,
                           QN("m",0)=>1,
                           QN("m",2)=>1,
                           QN("m",0)=>1,
                           QN("m",2)=>1,
                           QN("m",-2)=>1,
                           QN("m",0)=>1; tags="b$(i)"))
    end
    ρ[1] = ITensor(T00, ind_b[1], dag(ind_p[1]), ind_p[1]')
    for i in 2:N-1
        ρ[i] = ITensor(T01, dag(ind_b[i-1]), ind_b[i], dag(ind_p[i]), ind_p[i]')
    end
    ρ[N] = ITensor(T02, dag(ind_b[N-1]), dag(ind_p[N]), ind_p[N]')

    state = Int[]
    for i in 1:N
        if i % 2 == 0
            push!(state, 1)
        else
            push!(state, 2)
        end
    end
    psi0_init = randomMPS(ind_p,state,2)
    sweeps = Sweeps(200)
    setmaxdim!(sweeps, 10,10,20,20,50,50,100,100,100,100,
               200,200,200,200,300,300,300,400,400,400,
               500,500,500,500,800,800,800,800,800,1000)
    setnoise!(sweeps, 1E-7,1E-7,1E-7,1E-7,0)
    E0, psi0 = dmrg(-ρ, psi0_init, sweeps)
    
    SvN_l = zeros(N-1)

    orthogonalize!(psi0, 1)
    U,S,V = svd(psi0[1], (siteind(psi0,1)))
    SvN = 0.0
    for n in 1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)
    end
    SvN_l[1] = SvN

    for b in 2:N-1
        orthogonalize!(psi0, b)
        U,S,V = svd(psi0[b], (linkind(psi0, b-1), siteind(psi0,b)))
        SvN = 0.0
        for n in 1:dim(S, 1)
            p = S[n,n]^2
            SvN -= p * log(p)
        end
        SvN_l[b] = SvN
    end

    chi = sweeps.maxdim[end]
    open("./data/SvN_dmrg_mpo_S_$(SS/2)_L_$(N)_chi_$(chi).dat", "w") do io
        write(io, SvN_l)
    end

    orthogonalize!(psi0, 1)
    corr = zeros(N-1)

    for i in 2:N
        corr_i = op("sz",ind_p[1]) * psi0[1]
        corr_i = corr_i * dag(prime(psi0[1]))
        for j in 2:i-1
            corr_i = corr_i * psi0[j]
            corr_i = corr_i * dag(prime(psi0[j],"Link"))
        end
        corr_i = corr_i * psi0[i]
        corr_i = corr_i * op("sz",ind_p[i])
        corr_i = corr_i * dag(prime(prime(psi0[i],"p$(i)"),"Link,l=$(i-1)"))
        corr[i-1] = corr_i[1]
    end

    open("./data/corr_dmrg_mpo_S_$(SS/2)_L_$(N)_chi_$(chi).dat", "w") do io
        write(io, corr)
    end

    psi1_init = randomMPS(ind_p,state,2)
    E1, psi1 = dmrg(-ρ,[psi0],psi1_init, sweeps)

    psi2_init = randomMPS(ind_p,state,2)
    E2, psi2 = dmrg(-ρ,[psi0,psi1],psi2_init, sweeps)

    psi3_init = randomMPS(ind_p,state,2)
    E3, psi3 = dmrg(-ρ,[psi0,psi1,psi2],psi3_init, sweeps)
    #println([E0,E1,E2,E3])
    open("./data/E_dmrg_mpo_S_$(SS/2)_L_$(N)_chi_$(chi).dat", "w") do io
        write(io, [E0,E1,E2,E3])
    end
    return T0
end

N = parse(Int, ARGS[1])
SS = parse(Int, ARGS[2])
T0 = main(N, SS)
    
