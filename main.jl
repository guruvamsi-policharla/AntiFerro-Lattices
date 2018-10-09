#Includes
@everywhere include("skyrm_aux.jl")
@everywhere using Distributions
@everywhere using StatsBase

#driver
Tmin = 0.6
Tchange = 0.1
Tmax = 0.6
N = 4
Temperature = Tmin:Tchange:Tmax

J_space = 1:0.5:10

skyrm_temp = SharedArray{Float64,3}(length(Temperature),length(J_space),nprocs()-1)
skyrm_err_temp = SharedArray{Float64,3}(length(Temperature),length(J_space),nprocs()-1)

mag_temp = SharedArray{Float64,3}(length(Temperature),length(J_space),nprocs()-1)
mag_err_temp = SharedArray{Float64,3}(length(Temperature),length(J_space),nprocs()-1)

qFT = SharedArray{Float64,4}(N,N,length(J_space),nprocs()-1)
@parallel for i in 2:nprocs()
    skyrm_temp[:,:,i-1],skyrm_err_temp[:,:,i-1],mag_temp[:,:,i-1],mag_err_temp[:,:,i-1],qFT[:,:,:,i-1] = fetch(@spawnat i montecarlo(Temperature,N,J_space))
end
