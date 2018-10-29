#Includes
using Distributed
addprocs(Sys.CPU_THREADS)
#addprocs(8)
println(nprocs())
@everywhere using SharedArrays
@everywhere include("/home/guru/repos/antiFerro/skyrm_aux.jl")
#@everywhere include("skyrm_aux.jl")
@everywhere using Distributions
@everywhere using StatsBase
@everywhere using LinearAlgebra
@everywhere using JLD2
#driver
Tmin = 0.1
Tchange = 0.1
Tmax = 2
N = 4
Temperature = Tmin:Tchange:Tmax

J_space = [0.0:0.1:0.3;0.35:0.02:0.65;0.7:0.1:1]

skyrm_temp = SharedArray{Float64,5}(length(Temperature),length(J_space),4,3,nprocs()-1)
skyrm_err_temp = SharedArray{Float64,5}(length(Temperature),length(J_space),4,3,nprocs()-1)

mag_temp = SharedArray{Float64,5}(length(Temperature),length(J_space),4,3,nprocs()-1)
mag_err_temp = SharedArray{Float64,5}(length(Temperature),length(J_space),4,3,nprocs()-1)

qFT = SharedArray{Float64,4}(N,N,length(J_space),nprocs()-1)
proc_complete = SharedArray{Int,1}(nprocs())
for i in 1:nprocs()
    proc_complete[i] = 0
end
@distributed for i in 2:nprocs()
    skyrm_temp[:,:,:,:,i-1],skyrm_err_temp[:,:,:,:,i-1],mag_temp[:,:,:,:,i-1],mag_err_temp[:,:,:,:,i-1],qFT[:,:,:,i-1] = fetch(@spawnat i montecarlo(Temperature,N,J_space))
    proc_complete[i] = 1
    println(i)
end

proc_complete[1] = 1
for i in 1:5000
    if(mean(proc_complete) == 1)
        println(proc_complete)
        @save "data4x4fullres.jld2" skyrm_temp skyrm_err_temp mag_temp mag_err_temp Temperature N J_space
	break
    end
    println(proc_complete)
    sleep(350)
end

println("ending")
