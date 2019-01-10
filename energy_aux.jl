function total_mag(lat)
    #returns total magnetisation vector
    return sum(lat)
end

function test_flip(x, y, J, lat, T)
""" Checks whether energy allows for a flip or not """
    #a = sample_uni()
    a = sample_gauss(lat[x,y])
    de = -energy_pos(x,y,J,lat) + energy_pos(x,y,J,lat,a);

    if(de<0)
        lat[x,y] = a
        return true
    elseif(rand() < exp(-de/T))
        lat[x,y] = a
        return true
    else
        return false
    end
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

function total_energy(J,lat)
    e = 0.0
    for i = 1:size(lat,1)
        for j = 1:size(lat,2)
            e = e + energy_pos(i,j,J,lat)
        end
    end
    return e/2
end
