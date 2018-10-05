function sample_gauss(v)
    #Input is the central vector around which we flip
    var = 0.1
    v_new = [rand(Normal(v[1],var)),rand(Normal(v[2],var)),rand(Normal(v[3],var))]
    v_new = v_new/sqrt(vecdot(v_new,v_new))
    return convert(Vector,v_new)
end

function sample_uni()
    x = [rand(Uniform(-1,1)), rand(Uniform(-1,1)), rand(Uniform(-1,1))]
    x = x/sqrt(vecdot(x,x))
    return convert(Vector,x)
end

function initialise(M::Int, N::Int)
    """ Initialising the lattice with random values """
    #lat = Array{Float64, 3}(N, N)
    lat = Array{Vector{Float64},2}(M, N);
    for i = 1:M
        for j = 1:N
            lat[i,j] = sample_uni()
        end
    end
    return lat
end

function initialise_gauss(M::Int, N::Int)
    """ Initialising the lattice with random values """
    #lat = Array{Float64, 3}(N, N)
    lat = Array{Vector{Float64},2}(M, N);
    for i = 1:M
        for j = 1:N
            lat[i,j] = sample_gauss([1,0,0])
        end
    end
    return lat
end


function energy_pos(x, y, J, lat, a = [0,0,0])
    M = size(lat,1)
    N = size(lat,2)

    up = mod(y,N)+1
    down = mod(y-2,N)+1
    left = mod(x-2,M) + 1
    right = mod(x,M) + 1
    urc = [mod(x,M)+1,mod(y,N)+1]
    ulc = [mod(x-2,M)+1,mod(y,N)+1]
    lrc = [mod(x,M)+1,mod(y-2,N)+1]
    llc = [mod(x-2,M)+1,mod(y-2,N)+1]

    if(a == [0,0,0])
        energy = 1*vecdot(lat[x,y],(lat[left,y]+lat[right,y]+lat[x,up]+lat[x,down]));
        energy = energy + J*vecdot(lat[x,y],(lat[urc[1],urc[2]]+lat[ulc[1],ulc[2]]+lat[lrc[1],lrc[2]]+lat[llc[1],llc[2]]));
        return energy
    else
        energy = 1*vecdot(a,(lat[left,y]+lat[right,y]+lat[x,up]+lat[x,down]));
        energy = energy + J*vecdot(a,(lat[urc[1],urc[2]]+lat[ulc[1],ulc[2]]+lat[lrc[1],lrc[2]]+lat[llc[1],llc[2]]));
        return energy
    end

end

function test_flip(x, y, J, lat, T)
""" Checks whether energy allows for a flip or not """
    a = sample_uni()
    #a = sample_gauss(lat[x,y])
    de = -energy_pos(x,y,J,lat) + energy_pos(x,y,J,lat,a);

    if(de<0)
        lat[x,y] = a
        return true
    elseif(rand()<exp(-de/T))
        lat[x,y] = a
        return true
    else
        return false
    end
end

function transient_results(lat, transient::Int, J, T)
    """Takes lat as input and removes initial transients by running for transient number of steps"""
    M = size(lat,1)
    N = size(lat,2)
    for i = 1:transient
        for j = 1:M*N
                x = rand(1:M)
                y = rand(1:N)
                test_flip(x,y,J,lat,T)
        end
    end
end

function total_mag(lat)
    #returns total magnetisation vector
    return sum(lat)
end

function total_energy(J,lat)
    e = 0.0
    for i = 1:size(lat,1)
        for j = 1:size(lat,2)
            e = e + energy_pos(i,j,J,lat)
        end
    end
    return e/2
end

function exact_heisen_energy(T)
    #(e^(1/T) (-1 + T) T - e^(-1/T) T (1 + T))/(2 T sinh(1/T))
    return (exp(1/T) * (-1 + T) * T - exp(-1/T) * T * (1 + T))/(2*T*sinh(1/T))
end

function jackknife(vec)
    s = sum(vec)
    n = length(vec)
    vec_jack = (s - vec)/(n-1)
    jack_avg = sum(vec_jack) / n

    jack_err = sqrt(sum((vec_jack-jack_avg).^2) * (n-1)/n)
    #jack_err = sqrt((n-1)*(jack_err - jack_avg.^2))
    return jack_avg,jack_err
end

function four_trans(lat)
    x = zeros(M,N)
    y = zeros(M,N)
    z = zeros(M,N)
    for i in 1:M
      for j in 1:N
        x[i,j] = lat[i,j][1]
        y[i,j] = lat[i,j][2]
        z[i,j] = lat[i,j][3]
      end
    end

    ftx = abs.(fft(x))
    fty = abs.(fft(y))
    ftz = abs.(fft(z))

    ft = abs.(sqrt.(ftx.*ftx + fty.*fty + ftz.*ftz))
return ft
end

function skyrmion_number(lat)
    M = size(lat,1)
    N = size(lat,2)
    q = 0
    for i in 1:M
        for j in 1:N
            a = lat[i,j]#centre
            b = lat[i,mod(j,N)+1]#right
            c = lat[mod(i,M)+1,mod(j,N)+1]#rightdown
            d = lat[mod(i,M)+1,j]#down
            q = q + (spher_tri_area(a,b,c) + spher_tri_area(a,c,d))/(4*pi)
        end
    end
    return q
end

function spher_tri_area(a,b,c)
    x = cross(a,b)
    y = cross(b,c)
    z = cross(c,a)
    try
        a1 = acos(vecdot(x,-z)/vecnorm(x)/vecnorm(z))
        a2 = acos(vecdot(y,-x)/vecnorm(y)/vecnorm(x))
        a3 = acos(vecdot(z,-y)/vecnorm(z)/vecnorm(y))
        return (a1 + a2 + a3 - pi)*sign(vecdot(a,cross(b,c)))
    catch err
        if isa(err,DomainError)
            println(vecdot(x,-z)/vecnorm(x)/vecnorm(-z),vecdot(y,-x)/vecnorm(y)/vecnorm(-x),vecdot(z,-y)/vecnorm(z)/vecnorm(-y))
            a1 = acos(round(vecdot(x,-z)/vecnorm(x)/vecnorm(z)))
            a2 = acos(round(vecdot(y,-x)/vecnorm(y)/vecnorm(x)))
            a3 = acos(round(vecdot(z,-y)/vecnorm(z)/vecnorm(y)))
            return (a1 + a2 + a3 - pi)*sign(vecdot(a,cross(b,c)))
        end
    end
end


function montecarlo(Temperature,N,J_space)
    #Tmin = 0.0001
    #Tchange = 0.05
    #Tmax = 0.2
    mcs = 40000
    M = N

    norm=(1.0/float(M*N))

    #Temperature = Tmin:Tchange:Tmax

    M_vec = zeros(length(Temperature),2)
    JM_vec = zeros(length(Temperature),length(J_space))
    JM_vec_err = zeros(length(Temperature),length(J_space))
    M_jack = zeros(mcs,1)

    Jskyrm_vec = zeros(length(Temperature),length(J_space))
    Jskyrm_vec_err = zeros(length(Temperature),length(J_space))
    skyrm_vec = zeros(length(Temperature),2)
    skyrm_jack = zeros(mcs,1)

    autocor_vec = 0
    Jcount = 1
    for J in J_space
        lat = initialise(M,N)
        count = 1
        for T in Temperature
            transient_results(lat,3000,J,T)
            E = total_energy(J,lat)
            for i in 1:mcs
                for j in 1:M*N
                    x = rand(1:M)
                    y = rand(1:N)
                    E_0 = energy_pos(x,y,J,lat)
                    if(test_flip(x,y,J,lat,T))
                        E = E + energy_pos(x,y,J,lat) - E_0
                    end
                end
                skyrm_num = skyrmion_number(lat)
                skyrm_jack[i] = (skyrm_num).^2

                Mag = total_mag(lat)
                M_jack[i] = vecnorm(Mag)
            end
            skyrm_jack = skyrm_jack*norm
            skyrm_vec[count,1], skyrm_vec[count,2] = jackknife(skyrm_jack)

            M_jack = M_jack*norm
            M_vec[count,1], M_vec[count,2] = jackknife(M_jack)

            count = count + 1
            println(T)
            #if T == Tmin
            #    autocor_vec = autocor(skyrm_jack)
            #end
        end
    Jskyrm_vec[:,Jcount] = skyrm_vec[:,1]
    Jskyrm_vec_err[:,Jcount] = skyrm_vec[:,2]

    JM_vec[:,Jcount] = M_vec[:,1]
    JM_vec_err[:,Jcount] = M_vec[:,2]

    Jcount = Jcount + 1
    end

    return Jskyrm_vec,Jskyrm_vec_err,JM_vec,JM_vec_err
end
