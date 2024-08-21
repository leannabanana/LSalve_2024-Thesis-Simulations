include("methods/EI_estimate.jl")
include("methods/chaotic_system_methods.jl")

Random.seed!(1234)

n_orbits = 10^3
orbit_length = 10^6
initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)
window_sizes = collect(1:10)
a = 2   
pertubation = 1/10^3
p0 = 0
p1 = 1/(sqrt(2))



orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
average_orbits = observable_two.(orbits, p1)
mat_orb = reduce(hcat, average_orbits)
frequency_plots_min(mat_orb, window_sizes)


estimate = EI_window_min(mat_orb, window_sizes)

mat_orb



new_mu = leanna_mu_2.(window_sizes, estimate[1][1], estimate[2][1], estimate[4], estimate[4][1])

dfmatrix = DataFrame(mat_orb, :auto)
CSV.write("Data_csv/updated_data/trajectories_106_1sqrt2.csv", estimate)

df = DataFrame(location = estimate[1], scale = estimate[2], EI = estimate[3], shape = estimate[4])
CSV.write("Output_Images/18-08-2024/gevfit_min_1sqrt2_data.csv", df, delim=',', header=true)

g1 = scatter(window_sizes, estimate[1], xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2",  ms=3, ma=1)
pl.plot!(window_sizes, new_mu, label=L"μ_2 = μ_1 λ^{-(k - 1)}")


g2 = scatter(window_sizes, estimate[2], xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", ms=3, ma=1)
pl.plot!(window_sizes, leanna_mu_2(window_sizes, ), label=L"σ_2 = σ_1 λ^{\exp(2)k - 1}")
pl.plot!(window_sizes, estimate[2][1].* 2 .^(-1.0 .* ( window_sizes .- 1)), label=L"σ_2 = σ_1 λ^{-(k - 1)}")

g3 = scatter(window_sizes, estimate[3], xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
g4 = scatter(window_sizes, estimate[4], xticks=1:1:13, xlabel = " k ", ylabel = L"ξ", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)

gplot = pl.plot(g1, g2, layout = @layout([A B]), size=(700,400), plot_title="GEV Fit", plot_titlevspan=0.04)

savefig(gplot,"Output_Images/18-08-2024/mingevfit_lines_1sqrt2.pdf")


orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
average_orbits_0 = observable_two.(orbits, 0)


mat_orb = reduce(hcat, average_orbits)
gumbelestimate = gumbel_window_min(mat_orb, window_sizes)
e1 = scatter(window_sizes, gumbelestimate[1], xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
pl.plot!(window_sizes, estimate[1][1].*(1 ./ 2 .^(window_sizes.-1.0)))

e2 = scatter(window_sizes, gumbelestimate[2], xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
e3 = scatter(window_sizes, gumbelestimate[3], xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)

gumbelmin = pl.plot(e1, e2, e3, layout = @layout([A B ; C D]), size=(700,550), plot_title="Gumbel Fit", plot_titlevspan=0.04)
savefig(gumbelmin,"Output_Images/18-08-2024/gumbel_min_1sqrt2.pdf")

estimate_av = EI_window_av(mat_orb, window_sizes)

estimate[2][1]

model(k, a) = estimate_av[1][1] ./ k .+ (a[1] ./k).*log.(estimate_av[3] ./ estimate_av[3][1])
a0 = [1.0]
fit = curve_fit(model, window_sizes, estimate_av[1], a0)
param = fit.param


m1 = scatter(window_sizes, estimate_av[1], xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title=L"Moving av $x_0 = 0$", mc="tomato2", legend=false, ms=3, ma=1)
pl.plot!(window_sizes, estimate_av[1][1] ./ window_sizes .+ (38.2271057231807 ./window_sizes).*log.(estimate_av[3] ./ estimate_av[3][1]))




m2 = scatter(window_sizes, estimate_av[2], xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
pl.plot!(window_sizes, estimate_av[2][1] ./ window_sizes)

m3 = scatter(window_sizes, estimate_av[3], xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
m4 = scatter(window_sizes, estimate_av[4], xticks=1:1:13, xlabel = " k ", ylabel = L"ξ", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)

avfit = pl.plot(m1, m2, m3, m4, layout = @layout([A B ; C D]), size=(700,550), plot_title="GEV Fit", plot_titlevspan=0.04)
savefig(avfit,"Output_Images/updated_plots/gevfit_av__1sqrt2.pdf")

estimate[3]
gumbel_av = gumbel_window_av(mat_orb, window_sizes)


gm1 = scatter(window_sizes, gumbel_av[1], xticks=1:1:13, legend=true, label="simulated", xlabel = " k ", ylabel = L"\mu", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", ms=3, ma=1)
# pl.plot!(window_sizes, gumbel_av[1][1] ./ log.(2 .* (window_sizes)), label = L"μ_2 = \dfrac{μ_1}{\log(2k)}") 
# pl.plot!(window_sizes, gumbel_av[1][1] ./ log.(2 .^ (window_sizes)), label = L"μ_2 = \dfrac{μ_1}{\log(2^k)}") 
pl.plot!(window_sizes, leanna_mu_2.(window_sizes, gumbel_av[1][1], gumbel_av[2][1], gumbel_av[3], gumbel_av[3][1]), label = L"μ_2 = \dfrac{μ_1}{k}")




gm2 = scatter(window_sizes, gumbel_av[2], xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend = true, label = "simulated",  ms=3, ma=1)
pl.plot!(window_sizes, gumbel_av[2][1] ./ window_sizes, label= L"σ_2 = \dfrac{σ_1}{k}")

gm3 = scatter(window_sizes, gumbel_av[3], xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)

gavfit = pl.plot(gm1, gm2, layout = @layout([A B]), size=(700,400), plot_title="Gumbel Fit", plot_titlevspan=0.04)


savefig(gavfit,"Output_Images/18-08-2024/gumbelfit_fits_1sqrt2.pdf")

m4 = scatter(window_sizes, estimate_av[4], xticks=1:1:13, xlabel = " k ", ylabel = L"ξ", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
