Tmin = 0.1
Tchange = 0.1
Tmax = 3 #Change temp IN BOTH LOCATIONS!!!
Temperature = Tmin:Tchange:Tmax

J_space = [0,0.25,0.5,0.75,1.0,1.5,2.0]
skyrm = Array{Float64,3}(length(Temperature),length(J_space),2)
mag = Array{Float64,3}(length(Temperature),length(J_space),2)

mag[:,:,1] = reshape(mean(mag_temp,3),(size(mag_temp,1),size(mag_temp,2)))
mag[:,:,2] = reshape(sqrt.(sum(mag_err_temp.^2,3))./size(mag_err_temp,3),(size(mag_err_temp,1),size(mag_err_temp,2)))
skyrm[:,:,1] = reshape(mean(skyrm_temp,3),(size(skyrm_temp,1),size(skyrm_temp,2)))
skyrm[:,:,2] = reshape(sqrt.(sum(skyrm_err_temp.^2,3))./size(skyrm_err_temp,3),(size(skyrm_err_temp,1),size(skyrm_err_temp,2)))

figure()

for i in 1:length(J_space)
    #errorbar(Temperature,mean(data["mag"*string(i)],2),mean(data["mag_err"*string(i)],2))
    errorbar(Temperature,mag[:,i,1],yerr = mag[:,i,2],fmt="o",linestyle="-")
end
#ax = gca() # Get the handle of the current axis
title("Magnetisation Curve "*string(N)*"x"*string(N))
xlabel("Temperature")
ylabel("abs(mag)")
grid("on")
legend("J1/J2 = ".*string.(J_space))

figure()
for i in 1:length(J_space)
    errorbar(Temperature,skyrm[:,i,1],skyrm[:,i,2],fmt="o",linestyle="-")
    #errorbar(Temperature,skyrm[:,:,1],sqrt.(sum(data["skyrm_err"*string(i)].^2,2))./size(data["skyrm_err"*string(i)],2))
end
title("Skyrmion^2 Curve "*string(N)*"x"*string(N))
xlabel("Temperature")
ylabel("skyrm^2")
grid("on")
legend("J1/J2 = ".*string.(J_space))

figure()
for i in 1:5:length(Temperature)
    errorbar(J_space,skyrm[i,:,1],skyrm[i,:,2],fmt="o",linestyle="-")
end
title("Skyrmion^2 Curve "*string(N)*"x"*string(N))
xlabel("J1/J2")
ylabel("skyrm^2")
grid("on")
legend("T = ".*string.(Temperature[1:5:length(Temperature)]))
