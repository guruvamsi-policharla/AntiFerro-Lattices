function sample_gauss(v)
    #Input is the central vector around which we flip
    var = 0.2
    v_new = [rand(Normal(v[1],var)),rand(Normal(v[2],var)),rand(Normal(v[3],var))]
    v_new = v_new/sqrt(dot(v_new,v_new))
    return convert(Vector,v_new)
end

function sample_uni()
    x = [rand(Uniform(-1,1)), rand(Uniform(-1,1)), rand(Uniform(-1,1))]
    x = x/sqrt(dot(x,x))
    return convert(Vector,x)
end

function initialise(M::Int, N::Int)
    """ Initialising the lattice with random values """
    #lat = Array{Float64, 3}(N, N)
    lat = Array{Vector{Float64},2}(undef,M, N);
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
        energy = 1*dot(lat[x,y],(lat[left,y]+lat[right,y]+lat[x,up]+lat[x,down]));
        energy = energy + J*dot(lat[x,y],(lat[urc[1],urc[2]]+lat[ulc[1],ulc[2]]+lat[lrc[1],lrc[2]]+lat[llc[1],llc[2]]));
        return energy
    else
        energy = 1*dot(a,(lat[left,y]+lat[right,y]+lat[x,up]+lat[x,down]));
        energy = energy + J*dot(a,(lat[urc[1],urc[2]]+lat[ulc[1],ulc[2]]+lat[lrc[1],lrc[2]]+lat[llc[1],llc[2]]));
        return energy
    end

end

function test_flip(x, y, J, lat, T)
""" Checks whether energy allows for a flip or not """
    #a = sample_uni()
    a = sample_gauss(lat[x,y])
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

function jackknife(v)
    s = sum(v)
    n = length(v)
    vec_jack = (s .- v)/(n-1)
    jack_avg = mean(vec_jack)
    jack_err = sqrt((mean(vec_jack.^2) .- jack_avg.^2) * (n-1))
    return jack_avg,jack_err
end

function bindjack(vec4,vec2)
    n = length(vec4)
    s4 = sum(vec4)
    s2 = sum(vec2)
    vec_jack = ((s4 .- vec4)./(s2 .- vec2).^2) .* (n-1)

    jack_avg = mean(vec_jack)
    jack_err = sqrt(abs(mean(vec_jack.^2) .- jack_avg.^2) * (n-1))

    #println("Finished calculating Binder")
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

function four_trans_skyrm(lat)
    M = size(lat,1)
    N = size(lat,2)
    q = zeros(M,N)
    for i in 1:M
       for j in 1:N
           a = lat[i,j]#centre
           b = lat[i,mod(j,N)+1]#right
           c = lat[mod(i,M)+1,mod(j,N)+1]#rightdown
           d = lat[mod(i,M)+1,j]#down
           q[i,j] = (spher_tri_area(a,b,c) + spher_tri_area(a,c,d))/(4*pi)
       end
    end
    ft = abs.(fft(q))
return ft
end

function trans_skyrm(lat)
    M = size(lat,1)
    N = size(lat,2)
    q = zeros(M,N)
    for i in 1:M
       for j in 1:N
           a = lat[i,j]#centre
           b = lat[i,mod(j,N)+1]#right
           c = lat[mod(i,M)+1,mod(j,N)+1]#rightdown
           d = lat[mod(i,M)+1,j]#down
           q[i,j] = (spher_tri_area(a,b,c) + spher_tri_area(a,c,d))/(4*pi)
       end
    end
return q
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
        a1 = acos(dot(x,-z)/norm(x)/norm(z))
        a2 = acos(dot(y,-x)/norm(y)/norm(x))
        a3 = acos(dot(z,-y)/norm(z)/norm(y))
        return (a1 + a2 + a3 - pi)*sign(dot(a,cross(b,c)))
    catch err
        if isa(err,DomainError)
            println(dot(x,-z)/norm(x)/norm(-z),dot(y,-x)/norm(y)/norm(-x),dot(z,-y)/norm(z)/norm(-y))
            a1 = acos(round(dot(x,-z)/norm(x)/norm(z)))
            a2 = acos(round(dot(y,-x)/norm(y)/norm(x)))
            a3 = acos(round(dot(z,-y)/norm(z)/norm(y)))
            return (a1 + a2 + a3 - pi)*sign(dot(a,cross(b,c)))
        end
    end
end


function plotlat(lat,index1=0,index2=0)
    #fig = figure()
    w, h = figaspect(0.4)
    fig = figure(figsize=(w,h))
    M = size(lat,1)
    N = size(lat,2)
    X = repmat(1:M,M)
    Y = []
    for i in 1:M
        Y = vcat(Y,repmat([i],M))
    end

    Z = zeros(size(lat,2)^2)
    U = zeros(M*N)
    V = zeros(M*N)
    W = zeros(M*N)
    col = zeros(M*N)
    for i in 1:M
        for j in 1:N
            U[(i-1)*M+j] = lat[i,j][1]
            V[(i-1)*M+j] = lat[i,j][2]
            W[(i-1)*M+j] = lat[i,j][3]
            if(lat[i,j][3]>0)
                col[(i-1)*M+j] = 'g'
            else
                col[(i-1)*M+j] = 'b'
            end
        end
    end

    q = zeros(M,N)
    for i in 1:M
       for j in 1:N
           a = lat[i,j]#centre
           b = lat[i,mod(j,N)+1]#right
           c = lat[mod(i,M)+1,mod(j,N)+1]#rightdown
           d = lat[mod(i,M)+1,j]#down
           q[i,j] = (spher_tri_area(a,b,c) + spher_tri_area(a,c,d))/(4*pi)
       end
    end

    #contourf(X, Y, q)
    #colorbar()
    #subplot(121,projection="3d")
    #quiver(X,Y,Z,U,V,W)
    #zlim(-1,1)
    #title("Skyrmion Number = "*string.(skyrmion_number(lat)))
    subplot(121)
    contourf(q,vmin=-1,vmax=1)
    colorbar()
    title("Skyrmion Number = "* string.(skyrmion_number(lat)))
    subplot(122)
    quiver(X,Y,U,V,col)
    title("Skyrmion Number = "* string.(skyrmion_number(lat)))
    savefig("/home/vamsi/Github/quiveranim/" * string((index1-1)*M*M+index2) * ".png")
    close()
end

function lat_transform(lat,latindex)
    N = size(lat,1)
    if latindex == 2
        for ii in 1:N
            for jj in 1:N
                lat[ii,jj] = (-1)^(jj).*lat[ii,jj]
            end
        end
    elseif latindex == 3
        for ii in 1:N
            for jj in 1:N
                lat[ii,jj] = (-1)^(ii).*lat[ii,jj]
            end
        end
    elseif latindex == 4
        for ii in 1:N
            for jj in 1:N
                lat[ii,jj] = (-1)^(ii+jj).*lat[ii,jj]
            end
        end
    end
end

function montecarlo(Temperature,N,J_space)
    mcs = 50000
    M = N

    normalisation=(1.0/float(M*N))
    qFT = zeros(M,N,length(J_space))

    JM_vec = zeros(length(Temperature),length(J_space),4,3)
    JM_vec_err = zeros(length(Temperature),length(J_space),4,3)
    Jskyrm_vec = zeros(length(Temperature),length(J_space),4,3)
    Jskyrm_vec_err = zeros(length(Temperature),length(J_space),4,3)
    Jmagbind_vec = zeros(length(Temperature),length(J_space),4)
    Jskyrmbind_vec = zeros(length(Temperature),length(J_space),4)
    Jmagbind_vec_err = zeros(length(Temperature),length(J_space),4)
    Jskyrmbind_vec_err = zeros(length(Temperature),length(J_space),4)
#################################################################
    M_vec = zeros(length(Temperature),2,3,4)
    M_jack = zeros(mcs,3,4)

    skyrm_vec = zeros(length(Temperature),2,3,4)
    skyrm_jack = zeros(mcs,3,4)

    magbind_vec = zeros(length(Temperature),2,4)
    skyrmbind_vec = zeros(length(Temperature),2,4)

    #autocor_vec = 0
    Jcount = 1
    for J in J_space
        lat = initialise(M,N)
        Tcount = 1
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

                for latindex in 1:4
                    if latindex == 1
                        skyrm_num = skyrmion_number(lat)
                        Mag = total_mag(lat)
                    elseif latindex == 2
                        lat_transform(lat,latindex)
                        skyrm_num = skyrmion_number(lat)
                        Mag = total_mag(lat)
                        lat_transform(lat,latindex)
                    elseif latindex == 3
                        lat_transform(lat,latindex)
                        skyrm_num = skyrmion_number(lat)
                        Mag = total_mag(lat)
                        lat_transform(lat,latindex)
                    elseif latindex == 4
                        lat_transform(lat,latindex)
                        skyrm_num = skyrmion_number(lat)
                        Mag = total_mag(lat)
                        lat_transform(lat,latindex)
                    end

                    skyrm_jack[i,1,latindex] = abs(skyrm_num*normalisation)
                    skyrm_jack[i,2,latindex] = (skyrm_num*normalisation).^2
                    skyrm_jack[i,3,latindex] = (skyrm_num*normalisation).^4

                    M_jack[i,1,latindex] = (norm(Mag)*normalisation)
                    M_jack[i,2,latindex] = (norm(Mag)*normalisation).^2
                    M_jack[i,3,latindex] = (norm(Mag)*normalisation).^4

                end
                #qFT[:,:,Jcount] = qFT[:,:,Jcount] + four_trans_skyrm(lat)
            end

            for jj in 1:4
                for ii in 1:3
                    skyrm_vec[Tcount,1,ii,jj], skyrm_vec[Tcount,2,ii,jj] = jackknife(skyrm_jack[:,ii,jj])
                    M_vec[Tcount,1,ii,jj], M_vec[Tcount,2,ii,jj] = jackknife(M_jack[:,ii,jj])
                end
                magbind_vec[Tcount,1,jj],magbind_vec[Tcount,2,jj] = bindjack(M_jack[:,3,jj],M_jack[:,2,jj])
                skyrmbind_vec[Tcount,1,jj],skyrmbind_vec[Tcount,2,jj] = bindjack(skyrm_jack[:,3,jj],skyrm_jack[:,2,jj])
	    end
            Tcount = Tcount + 1

            #if T == Tmin
            #    autocor_vec = autocor(skyrm_jack)
            #end
        end

        for jj in 1:4
            for ii in 1:3
                Jskyrm_vec[:,Jcount,jj,ii] = skyrm_vec[:,1,ii,jj]
                Jskyrm_vec_err[:,Jcount,jj,ii] = skyrm_vec[:,2,ii,jj]

                JM_vec[:,Jcount,jj,ii] = M_vec[:,1,ii,jj]
                JM_vec_err[:,Jcount,jj,ii] = M_vec[:,2,ii,jj]
            end
            Jmagbind_vec[:,Jcount,jj] = magbind_vec[:,1,jj]
            Jmagbind_vec_err[:,Jcount,jj] = magbind_vec[:,2,jj]

            Jskyrmbind_vec[:,Jcount,jj] = skyrmbind_vec[:,1,jj]
            Jskyrmbind_vec_err[:,Jcount,jj] = skyrmbind_vec[:,2,jj]
        end

        Jcount = Jcount + 1
        println("J_space:",J)
    end

    return Jskyrm_vec,Jskyrm_vec_err,JM_vec,JM_vec_err,Jmagbind_vec,Jmagbind_vec_err,Jskyrmbind_vec,Jskyrmbind_vec_err,qFT
end
