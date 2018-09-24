#Includes
@everywhere include("skyrm_aux.jl")
@everywhere using Distributions
@everywhere using StatsBase
#driver
Tmin = 0.001
Tchange = 0.01
Tmax = 0.2 #Change temp IN BOTH LOCATIONS!!!
Temperature = Tmin:Tchange:Tmax

skyrm = SharedArray{Float64}(length(Temperature),nprocs()-1)
skyrm_err = SharedArray{Float64}(length(Temperature),nprocs()-1)
@parallel for i in 2:nprocs()
    skyrm[:,i-1],skyrm_err[:,i-1] = fetch(@spawnat i montecarlo(Temperature))
end
