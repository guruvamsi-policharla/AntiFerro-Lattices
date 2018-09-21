#Includes
include("aux.jl")
using Distributions
using Plots
using StatsBase

#Code
Tmin = 0.01
Tchange = 0.1
Tmax = 2.5
mcs = 100000
M = 4
N = 4



norm=(1.0/float(M*N))

Temperature = Tmin:Tchange:Tmax
J_space = [0]
E_vec = zeros(length(Temperature),2)
M_vec = zeros(length(Temperature),2)
JE_vec = zeros(length(Temperature),length(J_space))
JM_vec = zeros(length(Temperature),length(J_space))
Jskyrm_vec = zeros(length(Temperature),length(J_space))

skyrm_vec = zeros(length(Temperature),2)
acc_vec = zeros(length(Temperature),2)

E_jack = zeros(mcs,1)
M_jack = zeros(mcs,1)
acc_jack = zeros(mcs,1)
skyrm_jack = zeros(mcs,1)
autocor_vec = 0

gr()
for J in J_space
    lat = initialise(M,N)
    minFT = zeros(M,N)
    count = 1
    for T in Temperature
        transient_results(lat,3000,J,T)
        Mag = total_mag(lat)
        E = total_energy(J,lat)

        acc_rat = 0;
        for i in 1:mcs
            for j in 1:M*N
                x = rand(1:M)
                y = rand(1:N)
                E_0 = energy_pos(x,y,J,lat)
                Mag_0 = lat[x,y]
                skyrm_num_0 = skyrm_nn(x,y,lat)
                if(test_flip(x,y,J,lat,T))
                    acc_rat = acc_rat + 1
                    E = E + energy_pos(x,y,J,lat) - E_0
                    Mag = Mag + lat[x,y] - Mag_0
                end
            end
            skyrm_num = skyrmion_number(lat)
            mabs = sqrt(vecdot(Mag,Mag))
            E_jack[i] = E
            M_jack[i] = mabs
            acc_jack[i] = acc_rat
            skyrm_jack[i] = (skyrm_num)^2
        end
        E_jack = E_jack*norm
        M_jack = M_jack*norm
        acc_jack = acc_jack*norm
        skyrm_jack = skyrm_jack*norm
        M_vec[count,1], M_vec[count,2] = jackknife(M_jack)
        E_vec[count,1], E_vec[count,2] = jackknife(E_jack)
        acc_vec[count,1], acc_vec[count,2] = jackknife(acc_jack)
        skyrm_vec[count,1], skyrm_vec[count,2] = jackknife(skyrm_jack)
        println(T," ",skyrm_vec[count,1])
        count = count + 1
        if T == Tmin
            autocor_vec = autocor(skyrm_jack)
        end
        #if T == Tmin
        #    minFT = minFT + four_trans(lat)
        #end
    end
JM_vec[:,findfirst(J_space,J)] = M_vec[:,1]
JE_vec[:,findfirst(J_space,J)] = E_vec[:,1]
Jskyrm_vec[:,findfirst(J_space,J)] = skyrm_vec[:,1]
#l = 1:4:M
#xl = 2*(l-1)./M
#xl = string.(xl)
#xl = xl.*"pi"
#display(contourf(minFT,title = string(J),xticks = (l,xl),yticks = (l,xl)))
#savefig(string(J))
end
#jackknife analysis
