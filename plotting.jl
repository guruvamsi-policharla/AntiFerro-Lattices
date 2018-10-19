using PyPlot
using Statistics
using JLD2
#Tmin = 0.1
#Tchange = 0.1
#Tmax = 3 #Change temp IN BOTH LOCATIONS!!!
#Temperature = Tmin:Tchange:Tmax

#J_space = [0,0.25,0.5,0.75,1.0,1.5,2.0]
#J_space = 1:0.5:2.5

f=jldopen("./data/Full Fledged/data16x16full.jld2","r")
mag_temp = f["mag_temp"].s
skyrm_temp = f["skyrm_temp"].s
skyrm_err_temp = f["skyrm_err_temp"].s
mag_err_temp = f["mag_err_temp"].s
N = f["N"]
Temperature = f["Temperature"]
J_space = f["J_space"]

#skyrm = Array{Float64,3}(length(Temperature),length(J_space),2)
skyrm = Array{Float64,4}(undef,length(Temperature),length(J_space),4,3)
mag = Array{Float64,4}(undef,length(Temperature),length(J_space),4,3)
skyrm_err = Array{Float64,4}(undef,length(Temperature),length(J_space),4,3)
mag_err = Array{Float64,4}(undef,length(Temperature),length(J_space),4,3)

mag[:,:,:,:] = reshape(mean(mag_temp,dims=5),(size(mag_temp,1),size(mag_temp,2),4,3))
mag_err[:,:,:,:] = reshape(sqrt.(sum(mag_err_temp.^2,dims=5))./size(mag_err_temp,4),(size(mag_err_temp,1),size(mag_err_temp,2),4,3))
skyrm[:,:,:,:] = reshape(mean(skyrm_temp,dims=5),(size(skyrm_temp,1),size(skyrm_temp,2),4,3))
skyrm_err[:,:,:,:] = reshape(sqrt.(sum(skyrm_err_temp.^2,dims=5))./size(skyrm_err_temp,3),(size(skyrm_err_temp,1),size(skyrm_err_temp,2),4,3))


for jj in 1:3
figure()
    for ii in 1:4
        subplot(2,2,ii)
        for i in 1:length(J_space)
            errorbar(Temperature,mag[:,i,ii,jj],yerr = mag_err[:,i,ii,jj],fmt="o",linestyle="-")
            title(" Curve "*string(N)*"x"*string(N))
        end
        if ii == 1
            title("Magnetisation 00 - "*string(N)*"x"*string(N))
        elseif ii == 2
            title("Magnetisation 0pi - "*string(N)*"x"*string(N))
            legend("J1/J2 = ".*string.(J_space),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
        elseif ii == 3
            title("Magnetisation pi0 - "*string(N)*"x"*string(N))
        elseif ii == 4
            title("Magnetisation pipi - "*string(N)*"x"*string(N))
        end
        if ii>2
            xlabel("Temperature")
        end
        if mod(ii,2)==1
            if jj == 1
                ylabel("abs(mag)")
            elseif jj == 2
                ylabel("mag.^2")
            elseif jj == 3
                ylabel("mag.^4")
            end
        end
        grid("on")
    end
end

for jj in 1:3
    figure()
    for ii in 1:4
        subplot(2,2,ii)
        for i in 1:length(J_space)
            #errorbar(Temperature,mean(data["mag"*string(i)],2),mean(data["mag_err"*string(i)],2))
            errorbar(Temperature,skyrm[:,i,ii,jj],yerr = skyrm_err[:,i,ii,jj],fmt="o",linestyle="-")
        end
        if ii == 1
            title("Skyrmion 00 - "*string(N)*"x"*string(N))
        elseif ii == 2
            title("Skyrmion 0pi - "*string(N)*"x"*string(N))
            legend("J1/J2 = ".*string.(J_space),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
        elseif ii == 3
            title("Skyrmion pi0 - "*string(N)*"x"*string(N))
        elseif ii == 4
            title("Skyrmion pipi - "*string(N)*"x"*string(N))
        end
        if ii>2
            xlabel("Temperature")
        end
        if mod(ii,2)==1
            if jj == 1
                ylabel("abs(skyrm)")
            elseif jj == 2
                ylabel("skyrm.^2")
            elseif jj == 3
                ylabel("skyrm.^4")
            end
        end
        grid("on")
    end
end

for jj in 1:3
    figure()
    for ii in 1:4
        subplot(2,2,ii)
        for i in 1:2:length(Temperature)
            errorbar(J_space,skyrm[i,:,ii,jj],skyrm_err[i,:,ii,jj],fmt="o",linestyle="-")
        end
        if ii == 1
            title("Skyrmion 00 - "*string(N)*"x"*string(N))
        elseif ii == 2
            title("Skyrmion 0pi - "*string(N)*"x"*string(N))
            legend("T = ".*string.(Temperature[1:2:end]),bbox_to_anchor=[1.05,1],loc=2,ncol = 1)
        elseif ii == 3
            title("Skyrmion pi0 - "*string(N)*"x"*string(N))
        elseif ii == 4
            title("Skyrmion pipi - "*string(N)*"x"*string(N))
        end
        if ii>2
            xlabel("J1/J2")
        end
        if mod(ii,2)==1
            if jj == 1
                ylabel("abs(skyrm)")
            elseif jj == 2
                ylabel("skyrm.^2")
            elseif jj == 3
                ylabel("skyrm.^4")
            end
        end
        grid("on")
    end
end
