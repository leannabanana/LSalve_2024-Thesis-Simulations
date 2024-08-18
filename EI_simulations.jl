include("methods/chaotic_system_methods.jl")
include("methods/EI_estimate.jl")

Random.seed!(1234)

### Define initial_conditions of length n_orbits
n_orbits = 10^3
orbit_length = 10^5
initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)
window_sizes = collect(1:13)
a = 2   
pertubation = 1/10^3
p0 = 0
p1 = 1/(sqrt(2))

av_orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
observed_orbits_av = observable_two.(av_orbits, p1)
orbit_matrix_2 = reduce(hcat, observed_orbits_av)
estimate = EI_window_min(orbit_matrix_2, window_sizes)

df1 = DataFrame(location = estimate[1], scale = estimate[2], EI = estimate[3])
# CSV.write("Data_csv/min_nonrecurrent_windowsize_EI.csv", df1, delim=',', header=true)
g1 = scatter(window_sizes, estimate[1], xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
pl.plot!(window_sizes,  estimate[1][1] ./ ( exp(2) .^(window_sizes.-1.0)))


g2 = scatter(window_sizes, estimate[2], xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
g3 = scatter(window_sizes, estimate[3], xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)

savefig(g1,"Output_Images/22-07-2024/moving_minimum/mu_nonrecurrent.pdf")
savefig(g2,"Output_Images/22-07-2024/moving_minimum/sigma_nonrecurrent.pdf")
savefig(g3,"Output_Images/22-07-2024/moving_minimum/theta_nonrecurrent.pdf")

estimate_2 = EI_window_av(orbit_matrix_2, window_sizes)
df2 = DataFrame(location = estimate_2[1], scale = estimate_2[2], EI = estimate_2[3])
CSV.write("Data_csv/av_nonrecurrent_windowsize_EI.csv", df2, delim=',', header=true)

g4 = scatter(window_sizes, estimate_2[1], xlabel = " k ", ylabel = L"\mu", title=L"Moving average $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ma=1)
g5 = scatter(window_sizes, estimate_2[2], xlabel = " k ", ylabel = L"\sigma", title=L"Moving average $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false)
g6 = scatter(window_sizes, estimate_2[3], xlabel = " k ", ylabel = L"\theta", title=L"Moving average $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false)
savefig(g4,"Output_Images/22-07-2024/moving_average/mu_nonrecurrent.pdf")
savefig(g5,"Output_Images/22-07-2024/moving_average/sigma_nonrecurrent.pdf.pdf")
savefig(g6,"Output_Images/22-07-2024/moving_average/theta_nonrecurrent.pdf.pdf")


observed_orbits_av_1 = observable_two.(av_orbits, 0)
orbit_matrix_3 = reduce(hcat, observed_orbits_av_1)
estimate_5 = EI_window_min(orbit_matrix_3, window_sizes)
df3 = DataFrame(location = estimate_5[1], scale = estimate_5[2], EI = estimate_5[3])
CSV.write("Data_csv/min_recurrent_windowsize_EI.csv", df3, delim=',', header=true)

g7 = scatter(window_sizes, estimate_5[1], xlabel = " k ", ylabel = L"\mu", title=L"Moving min $x_0 = 0$", mc="tomato2", legend=false)
pl.plot!(window_sizes,  estimate_5[1][1] ./ ( exp(2) .^(window_sizes.-1.0)))

g8 = scatter(window_sizes, estimate_5[2], xlabel = " k ", ylabel = L"\sigma", title=L"Moving min $x_0 = 0$", mc="tomato2", legend=false)
g9 = scatter(window_sizes, estimate_5[3], xlabel = " k ", ylabel = L"\theta", title=L"Moving min $x_0 = 0$", mc="tomato2", legend=false)

savefig(g7,"Output_Images/22-07-2024/moving_minimum/mu_recurrent.pdf")
savefig(g8,"Output_Images/22-07-2024/moving_minimum/sigma_recurrent.pdf.pdf")
savefig(g9,"Output_Images/22-07-2024/moving_minimum/theta_recurrent.pdf.pdf")


estimate_6 = EI_window_av(orbit_matrix_3, window_sizes)
df4 = DataFrame(location = estimate_6[1], scale = estimate_6[2], EI = estimate_6[3])
CSV.write("Data_csv/av_recurrent_windowsize_EI.csv", df4, delim=',', header=true)

g10 = scatter(window_sizes, estimate_6[1], xlabel = " k ", ylabel = L"\mu", title=L"Moving av $x_0 = 0$", mc="tomato2", legend=false)
g11 = scatter(window_sizes, estimate_6[2], xlabel = " k ", ylabel = L"\sigma", title=L"Moving av $x_0 = 0$", mc="tomato2", legend=false)
g12 = scatter(window_sizes, estimate_6[3], xlabel = " k ", ylabel = L"\theta", title=L"Moving av $x_0 = 0$", mc="tomato2", legend=false)

savefig(g10,"Output_Images/22-07-2024/moving_average/mu_recurrent.pdf")
savefig(g11,"Output_Images/22-07-2024/moving_average/sigma_recurrent.pdf.pdf")
savefig(g12,"Output_Images/22-07-2024/moving_average/theta_recurrent.pdf.pdf")

av = moving_average_matrix(orbit_matrix_2, 10)
maxes = maximum(av, dims=1)[:]
av_orbits_2 = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
observed_orbits_av2 = observable_two.(av_orbits_2, 0)
orbit_matrix_3 = reduce(hcat, observed_orbits_av2)

av_par_2 = EI_estimation_average(orbit_matrix_3, 10)

x_axis_final =  collect(10:size(orbit_matrix_2)[1])

σ_av_orbit =scatter(x_axis_final, av_par[2], legend=false, xlabel = "Block Length", ylabel=L"\sigma^*", ms=1/2, ma =1/2, mc="indianred1", markerstrokecolor="indianred1", title =L"GEV Maxima sampled from $φ(T^n(x))$ with $x_0 = \dfrac{1}{\sqrt{2}}$
(Moving average)",  gridcolor=:gray19, gridalpha=1/2)
μ_av_orbit = scatter(x_axis_final, av_par[1], legend=false, xlabel = "Block Length", ylabel=L"\mu^*", ms=1/2, ma =1/2, mc="indianred1",  markerstrokecolor="indianred1",  title =L"GEV Maxima sampled from $φ(T^n(x))$ with$x_0 = \dfrac{1}{\sqrt{2}}$
(Moving average)", gridcolor=:gray19, gridalpha=1/2)
θ_av = scatter(x_axis_final, av_par[3], legend=false, xlabel = "Block Length", ylabel=L"θ", ms=1.5, ma =1/3, mc="indianred1",  markerstrokecolor="black",  title =L"GEV Maxima sampled from $φ(T^n(x))$ with $x_0 = \dfrac{1}{\sqrt{2}}$
(Moving average)", gridcolor=:gray19, gridalpha=1/2)

savefig(σ_av_orbit,"Output_Images/EI_estimation_orbit/moving_average/irrational_moving_av_sigma.pdf")
savefig(μ_av_orbit,"Output_Images/EI_estimation_orbit/moving_average/irrational_moving_av_mu.pdf")
savefig(θ_av,"Output_Images/EI_estimation_orbit/moving_average/irrational_moving_av_EI.pdf")


σ_av_orbit_1 =scatter(x_axis_final, av_par_2[2], legend=false, xlabel = "Block Length", ylabel=L"\sigma^*", ms=1/2, ma =1/2, mc="indianred1", markerstrokecolor="indianred1", title =L"GEV Maxima sampled from $φ(T^n(x))$ with $x_0 = 0$
(Moving average)",  gridcolor=:gray19, gridalpha=1/2)
μ_av_orbit_1 = scatter(x_axis_final, av_par_2[1], legend=false, xlabel = "Block Length", ylabel=L"\mu^*", ms=1/2, ma =1/2, mc="indianred1",  markerstrokecolor="indianred1",  title =L"GEV Maxima sampled from $φ(T^n(x))$ with $x_0 = 0$
(Moving average)", gridcolor=:gray19, gridalpha=1/2)
θ_av_1 = scatter(x_axis_final, av_par_2[3], legend=false, xlabel = "Block Length", ylabel=L"θ", ms=1.5, ma =1/3, mc="indianred1",  markerstrokecolor="black",  title =L"GEV Maxima sampled from $φ(T^n(x))$ with $x_0 = 0$
(Moving average)", gridcolor=:gray19, gridalpha=1/2)


savefig(σ_av_orbit_1,"Output_Images/EI_estimation_orbit/moving_average/periodic_moving_av_sigma.pdf")
savefig(μ_av_orbit_1,"Output_Images/EI_estimation_orbit/moving_average/periodic_moving_av_mu.pdf")
savefig(θ_av_1,"Output_Images/EI_estimation_orbit/moving_average/periodic_moving_av_EI.pdf")

