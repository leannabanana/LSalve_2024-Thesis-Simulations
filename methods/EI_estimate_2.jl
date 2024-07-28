include("chaotic_system_methods.jl")
include("EI_estimate.jl")

n_orbits = 10^3
orbit_length = 10^4
initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)

### Simulate orbits
a = 2 
pertubation = 1/10^3
p0 = 0
# p1 = 1/(sqrt(2*pi))


av_orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
observed_orbits_av = observable_two.(av_orbits, 1/sqrt(2))
orbit_matrix_2 = reduce(hcat, observed_orbits_av)
av_par = EI_estimation_average(orbit_matrix_2, 10)
data22 = DataFrame(location = av_par[1], scale =  av_par[2], EI =  av_par[3])
CSV.write("Data_csv/av_length_vs_parameters_nonrecurrent.csv", data22, delim=',', header=true)
av_par
av_orbits_5 = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
observed_orbits_av_5 = observable_two.(av_orbits_5, 0)
average_2 = reduce(hcat, observed_orbits_av_5)
av_par_2 = EI_estimation_average(average_2, 10)
data22 = DataFrame(location = av_par_2[1], scale =  av_par_2[2], EI =  av_par_2[3])
CSV.write("Data_csv/av_length_vs_parameters_recurrent.csv", data22, delim=',', header=true)



min_orbit = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
min_2 = observable_two.(min_orbit, 0)
min_mat_2 = reduce(hcat, min_2)
min_par = EI_estimation_min(min_mat_2, 10)
data33 = DataFrame(location = min_par[1], scale =  min_par[2], EI =  min_par[3])
CSV.write("Data_csv/min_length_vs_parameters_recurrent.csv", data33, delim=',', header=true)


min_orbit2 = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
minnn_2 = observable_two.(min_orbit2, 1/sqrt(2))
minn_2_matrix = reduce(hcat, minnn_2)
miiiii = EI_estimation_min(minn_2_matrix, 10)
data999 = DataFrame(location = miiiii[1][1:end-1], scale =  miiiii[2][1:end-1], EI =  miiiii[3][1:end-1])
CSV.write("Data_csv/min_length_vs_parameters_nonrecurrent.csv", data999, delim=',', header=true)

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

