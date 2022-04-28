#!/usr/bin/env julia

using PyPlot
using CorrelationFunctions
using Images
using FileIO
using Statistics
using Interpolations

im1 = load("../resize-effects/disks-small.png") .|> Bool
im2 = load("../resize-effects/disks-big.png")   .|> Bool
im3 = imresize(im1, (4000, 4000); method = BSpline(Linear()))

sv1 = Directional.surfvoid(im1, identity; periodic = true) |> mean
sv2 = Directional.surfvoid(im2, identity; periodic = true) |> mean
sv3 = Directional.surfvoid(im3, identity; periodic = true, void_phase = x -> (1 - x)^1.5) |> mean

figure(figsize = (10, 9), dpi = 300)
rc("font", size = 18)
plot(range(1, 2000, length(sv1)), sv1 * (300/4000))
plot(sv3)
plot(sv2)
legend([raw"$300\times 300$", raw"$4000 \times 4000$ (resize)", raw"$4000 \times 4000$ (calculated)"])
xlabel(raw"$r/a$")
ylabel(raw"$aF_{sv}(r)$")
ticklabel_format(scilimits = (0, 0))
savefig("../resize-effects/surfvoid.png")
