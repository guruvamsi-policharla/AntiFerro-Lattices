using PyPlot
using Statistics
using JLD2
#Tmin = 0.1
#Tchange = 0.1
#Tmax = 3 #Change temp IN BOTH LOCATIONS!!!
#Temperature = Tmin:Tchange:Tmax

#J_space = [0,0.25,0.5,0.75,1.0,1.5,2.0]
#J_space = 1:0.5:2.5
#=
f=jldopen("./data/data32bench1.jld2","r")
mag_temp = f["mag_temp"].s
skyrm_temp = f["skyrm_temp"].s
skyrm_err_temp = f["skyrm_err_temp"].s
mag_err_temp = f["mag_err_temp"].s
N = f["N"]
Temperature = f["Temperature"]
J_space = f["J_space"]
=#
#skyrm = Array{Float64,3}(length(Temperature),length(J_space),2)
skyrm = Array{Float64,3}(undef,length(Temperature),length(J_space),4)
mag = Array{Float64,3}(undef,length(Temperature),length(J_space),4)
skyrm_err = Array{Float64,3}(undef,length(Temperature),length(J_space),4)
mag_err = Array{Float64,3}(undef,length(Temperature),length(J_space),4)

mag[:,:,:] = reshape(mean(mag_temp,dims=4),(size(mag_temp,1),size(mag_temp,2),4))
mag_err[:,:,:] = reshape(sqrt.(sum(mag_err_temp.^2,dims=4))./size(mag_err_temp,4),(size(mag_err_temp,1),size(mag_err_temp,2),4))
skyrm[:,:,:] = reshape(mean(skyrm_temp,dims=4),(size(skyrm_temp,1),size(skyrm_temp,2),4))
skyrm_err[:,:,:] = reshape(sqrt.(sum(skyrm_err_temp.^2,dims=4))./size(skyrm_err_temp,3),(size(skyrm_err_temp,1),size(skyrm_err_temp,2),4))

figure()
for ii in 1:4
    subplot(2,2,ii)
    for i in 1:length(J_space)
        #errorbar(Temperature,mean(data["mag"*string(i)],2),mean(data["mag_err"*string(i)],2))
        errorbar(Temperature,mag[:,i,ii],yerr = mag_err[:,i,ii],fmt="o",linestyle="-")
        title("Magnetisation Curve "*string(N)*"x"*string(N))
        xlabel("Temperature")
        ylabel("abs(mag)")
        grid("on")
        legend("J1/J2 = ".*string.(J_space))
    end
end
#ax = gca() # Get the handle of the current axis


figure()
for ii in 1:4
    subplot(2,2,ii)
    for i in 1:length(J_space)
        #errorbar(Temperature,mean(data["mag"*string(i)],2),mean(data["mag_err"*string(i)],2))
        errorbar(Temperature,skyrm[:,i,ii],yerr = skyrm_err[:,i,ii],fmt="o",linestyle="-")
        title("Skyrmion^2 Curve "*string(N)*"x"*string(N))
        xlabel("Temperature")
        ylabel("skyrm^2")
        grid("on")
        legend("J1/J2 = ".*string.(J_space))
    end
end

figure()
for ii in 1:4
    subplot(2,2,ii)
    for i in 1:2:length(Temperature)
        errorbar(J_space,skyrm[i,:,ii],skyrm_err[i,:,ii],fmt="o",linestyle="-")
        title("Skyrmion^2 Curve "*string(N)*"x"*string(N))
        xlabel("J1/J2")
        ylabel("skyrm^2")
        grid("on")
        legend("T = ".*string.(Temperature[1:2:end]))
    end
end
