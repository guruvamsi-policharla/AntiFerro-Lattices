function total_mag(lat)
    #returns total magnetisation vector
    return sum(lat)
end

function myplus(a,b,N)
    c = mod(a+b,N)
    if c == 0
        c = N
    end
    return c
end

function test_flip(x, y, J, lat, T)
""" Checks whether energy allows for a flip or not """
    a0 = lat[x,y]
    a = sample_gauss(lat[x,y])
    de = -energy_pos(x,y,J,lat)
    lat[x,y] = a
    de = de + energy_pos(x,y,J,lat);

    if(de<0)
        lat[x,y] = a
        return true,de
    elseif(rand() < exp(-de/T))
        lat[x,y] = a
        return true
    else
        lat[x,y] = a0
        return false
    end
end

function energy_pos(x, y, J, lat, a = [0,0,0])
    M = size(lat,1)
    N = size(lat,2)

    up = lat[myplus(x,-1,N),y]
    down = lat[myplus(x,1,N),y]
    left = lat[x,myplus(y,-1,N)]
    right = lat[x,myplus(y,1,N)]
    urc = lat[myplus(x,-1,N),myplus(y,1,N)]
    ulc = lat[myplus(x,-1,N),myplus(y,-1,N)]
    lrc = lat[myplus(x,1,N),myplus(y,1,N)]
    llc = lat[myplus(x,1,N),myplus(y,-1,N)]

    if(a == [0,0,0])
        energy = 1*dot(lat[x,y], left+right+up+down );
        energy = energy + J*dot(lat[x,y],urc+ulc+lrc+llc);
        return energy
    else
        energy = 1*dot(a, left+right+up+down );
        energy = energy + J*dot(a,urc+ulc+lrc+llc);
        return energy
    end
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
