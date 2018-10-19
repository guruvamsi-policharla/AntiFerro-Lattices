#Includes
using Distributed
addprocs(Sys.CPU_THREADS)
println(nprocs())
@everywhere using SharedArrays
#@everywhere include("/home/guru/repos/antiFerro/skyrm_aux.jl")
@everywhere include("skyrm_aux.jl")
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

J_space = [0,0.25,0.5,0.75,1.0,1.5,2.0]

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
end

proc_complete[1] = 1
for i in 1:2000
    if(mean(proc_complete) == 1)
        println(proc_complete)
        @save "data64x64full.jld2" skyrm_temp skyrm_err_temp mag_temp mag_err_temp Temperature N J_space
	break
    end
    println(proc_complete)
    sleep(250)
end

println("ending")
