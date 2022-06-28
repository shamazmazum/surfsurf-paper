gencenters1(n) = [rand(Float64, 2) for _ in 1:n]

function gendisks1(side, centers, radius)
    array = zeros(Bool, (side, side))
    indices = CartesianIndices(array)

    for center in centers
        for idx in indices
            x = (idx[1] - 1) / side
            y = (idx[2] - 1) / side
            if (x - center[1])^2 + (y - center[2])^2 <= radius^2
                array[idx] = true
            end
        end
    end

    return array
end

function gencenters(side, λ)
    n = (λ * side^2) |> Poisson |> rand
    return reduce(hcat, (rand(1:side, 2) for i in 1:n))
end

function gendisks(side, R, λ)
    spheres = zeros(Bool, (side + 2R + 1, side + 2R + 1))
    sphere  = zeros(Bool, (2R + 1, 2R + 1))
    centers = gencenters(side, λ)
    for i in -R:R
        for j in -R:R
            dist = i^2 + j^2
            if dist < R^2
                sphere[j+R+1, i+R+1] = 1
            end
        end
    end
    for center in (centers[:,i] for i in 1:size(centers,2))
        x = center[1]
        y = center[2]
        spheres[x:x + 2R, y:y + 2R] .|= sphere
    end
    return spheres[R+1:end-R-1, R+1:end-R-1]
end

function ss_theory(r, R, λ)
    if r < 2R
        A = 4R^2 - r^2
        B = acos(r/(2R))
        part = exp(-λ*(r*sqrt(A)/2 + 2π*R^2 - 2*B*R^2))
        return part * ((4A*B^2 - 8π*A*B + 4*π^2*A)*R^2*λ^2*r + 4*sqrt(A)*R^2*λ) / (A*r)
    else
        return (2π*λ*R)^2*exp(-2π*λ*R^2)
    end
end

function sv_theory(r, R, λ)
     if r < 2R
        A = 4R^2 - r^2
        B = acos(r/(2R))
        return 2R*λ*(π-B)exp(-sqrt(A)*r*λ/2 + 2B*R^2*λ - 2π*R^2*λ)
    else
        return 2π*R*λ*exp(-2π*λ*R^2)
    end
end
