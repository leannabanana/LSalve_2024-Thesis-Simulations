include("methods/EI_estimate.jl")
include("methods/chaotic_system_methods.jl")

Random.seed!(1234)
n_orbits = 10^3
orbit_length = 10^5
initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)
window_sizes = collect(1:10)
a = 2   
pertubation = 1/10^3
p0 = 0
p1 = 1/(sqrt(2))

orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
average_orbits = observable_two.(orbits, 0)
mat_orb = reduce(hcat, average_orbits)
estimate = EI_window_min(mat_orb, window_sizes)

dfmatrix = DataFrame(mat_orb, :auto)
CSV.write("Data_csv/updated_data/trajectories_1sqrt2.csv", dfmatrix)

df = DataFrame(location = estimate[1], scale = estimate[2], EI = estimate[3], shape = estimate[4])
CSV.write("Output_Images/18-08-2024/gevfit_min_1sqrt2_data.csv", df, delim=',', header=true)
g1 = scatter(window_sizes, estimate[1], xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
pl.plot!(window_sizes, 0.5*(estimate[1][1]./ window_sizes))
g2 = scatter(window_sizes, estimate[2], xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
g3 = scatter(window_sizes, estimate[3], xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
g4 = scatter(window_sizes, estimate[4], xticks=1:1:13, xlabel = " k ", ylabel = L"ξ", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)

gplot = pl.plot(g1, g2, g3, g4, layout = @layout([A B ; C D]), size=(700,550), plot_title="GEV Fit", plot_titlevspan=0.04)
savefig(gplot,"Output_Images/18-08-2024/gevfit_min_1sqrt2.pdf")


orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
average_orbits = observable_two.(orbits, 0)
mat_orb = reduce(hcat, average_orbits)
gumbelestimate = gumbel_window_min(mat_orb, window_sizes)
e1 = scatter(window_sizes, gumbelestimate[1], xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
pl.plot!(window_sizes, estimate[1][1].*(1 ./ 2 .^(window_sizes.-1.0)))

e2 = scatter(window_sizes, gumbelestimate[2], xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
e3 = scatter(window_sizes, gumbelestimate[3], xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)

gumbelmin = pl.plot(e1, e2, e3, layout = @layout([A B ; C D]), size=(700,550), plot_title="Gumbel Fit", plot_titlevspan=0.04)
savefig(gumbelmin,"Output_Images/18-08-2024/gumbel_min_1sqrt2.pdf")

estimate_av = EI_window_av(mat_orb, window_sizes)
m1 = scatter(window_sizes, estimate_av[1], xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title=L"Moving av $x_0 = 0$", mc="tomato2", legend=false, ms=3, ma=1)
pl.plot!(window_sizes, (estimate_av[1][1] ./ log.(2 .* (window_sizes .- 1))))

m2 = scatter(window_sizes, estimate_av[2], xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
m3 = scatter(window_sizes, estimate_av[3], xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
m4 = scatter(window_sizes, estimate_av[4], xticks=1:1:13, xlabel = " k ", ylabel = L"ξ", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)

avfit = pl.plot(m1, m2, m3, m4, layout = @layout([A B ; C D]), size=(700,550), plot_title="GEV Fit", plot_titlevspan=0.04)
savefig(avfit,"Output_Images/18-08-2024/gevfit_av__1sqrt2.pdf")


gumbel_av = gumbel_window_av(mat_orb, window_sizes)
gm1 = scatter(window_sizes, gumbel_av[1], xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
gm2 = scatter(window_sizes, gumbel_av[2], xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
gm3 = scatter(window_sizes, gumbel_av[3], xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)

gavfit = pl.plot(gm1, gm2, gm3, layout = @layout([A B ; C D]), size=(700,550), plot_title="Gumbel Fit", plot_titlevspan=0.04)
savefig(gavfit,"Output_Images/18-08-2024/gumbelfit_av_1sqrt2.pdf")

m4 = scatter(window_sizes, estimate_av[4], xticks=1:1:13, xlabel = " k ", ylabel = L"ξ", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
