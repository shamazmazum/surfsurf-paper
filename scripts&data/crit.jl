#!/usr/bin/env julia

using PyPlot
using CorrelationFunctions

gencenters(n) = rand(Float64, (2, n))

function gendisks(centers, n, radius)
    disks = zeros(Bool, (n, n))
    for k in 1:size(centers, 2)
        center = centers[:, k]
        for i in 1:n
            for j in 1:n
                x = (i - 1) / n
                y = (j - 1) / n
                if (x - center[1])^2 + (y - center[2])^2 < radius^2
                    disks[j, i] = true
                end
            end
        end
    end
    return disks
end

centers = gencenters(700)
x = 100:100:4000
c = [gendisks(centers, n, 0.02) |> lowfreq_energy_ratio for n in x]

figure(figsize = (10, 9), dpi = 300)
rc("font", size = 18)
plot(collect(x), c, "b.")
xlabel(raw"$s$")
ylabel(raw"$C_{0.5}(s)$")
savefig("../images/plot-criterion.png")
