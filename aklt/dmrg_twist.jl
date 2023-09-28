using ITensors

function main(N::Int, SS::Int)
    ρ = MPO(N)
    T0 = reshape(collect(Float64,readeach(open("./data/mpo_M_aklt_S_$(SS/2)_L_$(N).dat", "r"), Float64)), 2,2,4,4) 
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
    sweeps = Sweeps(100)
    setmaxdim!(sweeps, 10,10,20,20,50,50,100,100,100,100,
               200,200,200,200,300,300,300)#,400)#,400,400,
               #500)#,500,500,500,800)
    setnoise!(sweeps, 1E-7,1E-7,1E-7,1E-7,0)
    E0, psi0 = dmrg(-ρ, psi0_init, sweeps)
    
    for i in 1:N
        orthogonalize!(psi0, i)
        u_i = ITensor([exp(1im*π*i/N) 0; 0 exp(-1im*π*i/N)], dag(ind_p[i]), ind_p[i]')
        p_i = u_i*psi0[i]
        noprime!(p_i)
        psi0[i] = p_i
    end
    psi1 = apply(ρ,psi0)
    p1 = inner(psi0,psi1)
    p2 = inner(psi1,psi1)
    println([-E0,p1,p2])
    chi = sweeps.maxdim[end]
    open("./data/twist_dmrg_mpo_S_$(SS/2)_L_$(N)_chi_$(chi)_rho.dat", "w") do io
        write(io, [-E0,p1,p2])
    end
    return T0
end

N = parse(Int, ARGS[1])
SS = parse(Int, ARGS[2])
T0 = main(N, SS)
    
