#!/usr/bin/env julia

using PyPlot
using CorrelationFunctions
using Images
using FileIO
using Statistics
using Interpolations

im1 = load("../resize-effects/disks-small.png") .|> Bool
im2 = load("../resize-effects/disks-big.png") .|> Bool
im3 = imresize(im1, (4000, 4000); method = BSpline(Linear()))

sv1 = Directional.surfsurf(im1, identity; periodic = true) |> mean
sv2 = Directional.surfsurf(im2, identity; periodic = true) |> mean
sv3 = Directional.surfsurf(im3, identity; periodic = true) |> mean

figure(figsize = (10, 9), dpi = 300)
rc("font", size = 18)
plot(range(1, 2000, length(sv1)), sv1 * (300/4000)^2)
plot(sv3)
plot(sv2)
xlim([100, 2000])
ylim([0, 0.0001])
legend([raw"$300\times 300$", raw"$4000 \times 4000$ (resize)", raw"$4000 \times 4000$ (calculated)"])
xlabel(raw"$r/a$")
ylabel(raw"$a^2F_{ss}(r)$")
ticklabel_format(scilimits = (0, 0))
savefig("../resize-effects/surfsurf.png")
