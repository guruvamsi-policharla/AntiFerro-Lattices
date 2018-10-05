#Includes
@everywhere include("skyrm_aux.jl")
@everywhere using Distributions
@everywhere using StatsBase

#driver
Tmin = 0.1
Tchange = 0.1
Tmax = 3
N = 8
Temperature = Tmin:Tchange:Tmax

J_space = [0,0.25,0.5,0.75,1.0,1.5,2.0]

skyrm_temp = SharedArray{Float64,3}(length(Temperature),length(J_space),nprocs()-1)
skyrm_err_temp = SharedArray{Float64,3}(length(Temperature),length(J_space),nprocs()-1)

mag_temp = SharedArray{Float64,3}(length(Temperature),length(J_space),nprocs()-1)
mag_err_temp = SharedArray{Float64,3}(length(Temperature),length(J_space),nprocs()-1)

@parallel for i in 2:nprocs()
    skyrm_temp[:,:,i-1],skyrm_err_temp[:,:,i-1],mag_temp[:,:,i-1],mag_err_temp[:,:,i-1] = fetch(@spawnat i montecarlo(Temperature,N,J_space))
end
