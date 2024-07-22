include("methods/chaotic_system_methods.jl")
include("methods/EI_estimate.jl")

Random.seed!(1234)

### Define initial_conditions of length n_orbits
n_orbits = 100
orbit_length = 100
initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)
window_sizes = collect(1:15)
### Simulate orbits
a = 2   
pertubation = 1/10^3
p0 = 0
p1 = 1/(sqrt(2))

av_orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
observed_orbits_av = observable_two.(av_orbits, 1/sqrt(2))
orbit_matrix_2 = reduce(hcat, observed_orbits_av)

av = moving_average_matrix(orbit_matrix_2, 10)
maxes = maximum(av, dims=1)[:]
extremal_FerroSegers(maxes, 0.95)

windowey = EI_window_av(orbit_matrix_2, window_sizes)
pl.plot(window_sizes, windowey[3])
av_par = EI_estimation_average(orbit_matrix_2, 10)


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

