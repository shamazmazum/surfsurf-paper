module DMvsSobel
using Distributions
using FFTW
using Images
using OffsetArrays
using Statistics
using Random
using PyPlot

# Disk generation code is copied from CorrelationFunctions.jl test suite
export gendisks, ss_theory, sv_theory,
    gradient, distance_map_edge, average,
    produce_fucking_graphics!

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

autocorr(arr :: AbstractArray) =
    (arr |> fft .|> abs2 |> ifft |> real) / length(arr)

function array_with_zero_based_indices(array :: Array)
    ax = map(x -> x .- 1, axes(array))
    return OffsetArray(array, ax)
end

function reflect(arr :: AbstractArray{T}) where T <: Real
    ft = arr |> fft |> array_with_zero_based_indices

    for idx in CartesianIndices(ft)
        ft[idx] = (-1)^(idx |> Tuple |> sum) * ft[idx]
    end

    return ft.parent |> ifft |> real
end

function gradient(arr :: AbstractMatrix, kernel = KernelFactors.sobel)
    x, y = imgradients(arr, kernel)
    return sqrt.(x.^2 + y.^2)
end

function distance_map_edge(arr :: AbstractMatrix{Bool})
    d = arr |> feature_transform |> distance_transform
    return @. 0 < d <= 1.42
end

function average(arr :: AbstractArray)
    arr = array_with_zero_based_indices(arr)
    center  = CartesianIndex(size(arr) .÷ 2)
    indices = CartesianIndices(arr)
    dict = Dict{Float64, Array{Float64}}()

    for idx in indices
        dist = sqrt(sum(Tuple(idx - center) .^2))
        vals = get(dict, dist, Float64[])
        push!(vals, arr[idx])
        dict[dist] = vals
    end

    xs = dict |> keys |> collect |> sort
    ys = map(xs) do key mean(dict[key]) end

    return xs, ys
end

function produce_fucking_graphics!()
    Random.seed!(1)
    img = DMvsSobel.gendisks(1000, 10, 0.002)
    save("surfsurf-paper/images/original.png", img[1:100, 1:100])

    function plot_it!(fn)
        interface = fn(img)
        cf = interface |> autocorr |> reflect
        xs, ys = average(cf)
        plot(xs, ys, linewidth = 2.0)
        return xs, cf
    end

    figure(figsize = (10, 9), dpi = 300)
    rc("font", size = 18)
    xs, cfs = plot_it!(gradient)
    xs, cfd = plot_it!(distance_map_edge)
    th(x) = ss_theory(x, 10, 0.002)
    plot(xs, th.(xs), linewidth = 2.0)
    xlabel(raw"$r$")
    ylabel(raw"$F_{ss}(r)$")
    legend(["Sobel kernel", "Distance transform", "Theory"]; loc = 2)
    xlim([0, 40])
    savefig("surfsurf-paper/images/Fss_mean_dir.png")
    save("surfsurf-paper/images/distance_map.png", img[1:100, 1:100] |> distance_map_edge)
    save("surfsurf-paper/images/sobel.png", img[1:100, 1:100] |> gradient)

    function plot_it!(cf, name)
        figure(figsize = (6, 6), dpi = 300)
        axis("off")
        s = CartesianIndex(size(cf) .÷ 2)
        w = CartesianIndex(50, 50)
        imshow(cf[(s - w):(s + w)])
        savefig("surfsurf-paper/images/Fss_map_$(name).png";
                bbox_inches = "tight",
                pad_inches = 0)
    end

    plot_it!(cfs, "sobel")
    plot_it!(cfd, "dm")
end

end
