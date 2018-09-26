#Includes
@everywhere include("skyrm_aux.jl")
@everywhere using Distributions
@everywhere using StatsBase
#driver
Tmin = 0.01
Tchange = 0.1
Tmax = 3 #Change temp IN BOTH LOCATIONS!!!
Temperature = Tmin:Tchange:Tmax

skyrm = SharedArray{Float64}(length(Temperature),nprocs()-1)
skyrm_err = SharedArray{Float64}(length(Temperature),nprocs()-1)


mag = SharedArray{Float64}(length(Temperature),nprocs()-1)
mag_err = SharedArray{Float64}(length(Temperature),nprocs()-1)
@parallel for i in 2:nprocs()
    skyrm[:,i-1],skyrm_err[:,i-1],mag[:,i-1],mag_err[:,i-1] = fetch(@spawnat i montecarlo(Temperature))
end
