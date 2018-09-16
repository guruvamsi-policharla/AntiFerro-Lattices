#Includes
include("aux.jl")
using Distributions
using Plots

#Code
Tmin = 0.01
Tchange = 0.1
Tmax = 4.0
mcs = 50000
M = 7
N = 7


lat = initialise(M,N)
minFT = zeros(M,N)
norm=(1.0/float(M*N))

Temperature = Tmin:Tchange:Tmax

E_vec = zeros(length(Temperature),2)
M_vec = zeros(length(Temperature),2)
JE_vec = zeros(length(Temperature),4)
JM_vec = zeros(length(Temperature),4)

acc_vec = zeros(length(Temperature),2)

E_jack = zeros(mcs,1)
M_jack = zeros(mcs,1)
acc_jack = zeros(mcs,1)

J_space = [0.0]

for J in J_space
    count = 1
    for T in Temperature

        transient_results(lat,1000,J,T)
        Mag = total_mag(lat)
        E = total_energy(J,lat)

        etot=0;
        mabstot=0;
        acc_rat = 0;

        for i in 1:mcs
            for j in 1:M*N
                x = rand(1:M)
                y = rand(1:N)
                E_0 = energy_pos(x,y,J,lat)
                Mag_0 = lat[x,y]
                if(test_flip(x,y,J,lat,T))
                    acc_rat = acc_rat + 1
                    E = E + energy_pos(x,y,J,lat) - E_0
                    Mag = Mag + lat[x,y] - Mag_0
                end
            end
            mabs = sqrt(vecdot(Mag,Mag))

            E_jack[i] = E
            M_jack[i] = mabs
            acc_jack[i] = acc_rat

            etot= etot + E
            mabstot = mabstot + mabs
        end
        E_jack = E_jack*norm
        M_jack = M_jack*norm
        acc_jack = acc_jack*norm

        M_vec[count,1], M_vec[count,2] = jackknife(M_jack)
        E_vec[count,1], E_vec[count,2] = jackknife(E_jack)
        acc_vec[count,1], acc_vec[count,2] = jackknife(acc_jack)
        println(T," ",M_vec[count,1])
        count = count + 1

        if T == Tmin
            minFT = minFT + four_trans(lat)
        end
    end
JM_vec[:,findfirst(J_space,J)] = M_vec[:,1]
JE_vec[:,findfirst(J_space,J)] = E_vec[:,1]
end
#jackknife analysis
