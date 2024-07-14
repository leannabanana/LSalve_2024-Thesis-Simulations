"""
This file plots parameters vs block length for our orbits
"""

include("methods/chaotic_system_methods.jl")
include("methods/parameter_methods.jl")

Random.seed!(1234)

### Define the window sizes
window_sizes = collect(1:50)

# Variables
orbit_length = 10^3
n_orbits = 10^3
block_length = orbit_length ./ window_sizes
initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)
x0 = 0
α = 1/3

# Parameters
a = 2 
perturbation = 1/10^3

# Simulate orbits
orbits = simulate_orbits(initial_conditions, a, orbit_length, perturbation, n_orbits)

param_block_1 = frechet_params_min(orbits, window_sizes, x0, α)

xis = pl.plot(block_length, param_block_1[1], title=L"ξ", legend=false)
mus = pl.plot(block_length, param_block_1[2], title=L"μ", legend=false)
theta = pl.plot(block_length, param_block_1[3], title=L"σ", legend=false)
min_params_1 = pl.plot(xis, mus, theta, size=(1000, 500), layout=(1,3), plot_title="Moving Minimum Parameters vs Window Size (Frechet)")

param_block_1 = frechet_params_av(orbits, window_sizes, x0, α)


xis = pl.plot(block_length, param_block_1[1], title=L"ξ", legend=false)
mus = pl.plot(block_length, param_block_1[2], title=L"μ", legend=false)
theta = pl.plot(block_length, param_block_1[3], title=L"σ", legend=false)
min_params_1 = pl.plot(xis, mus, theta, size=(1000, 500), layout=(1,3), plot_title="Moving Minimum Parameters vs Window Size (Frechet)")

parm_block_2 = gumbel_params_min(orbits, window_sizes, x0)

mus_2 = scatter(block_length, parm_block_2[2], title=L"μ", legend=false)
pl.plot!(block_length, parm_block_2[2][1] .+ parm_block_2[3][1].*log.(block_length))
theta_2 = pl.plot(block_length, parm_block_2[3], title=L"θ", legend=false)
min_params_2 = pl.plot(mus_2, theta_2, size=(1000, 500), layout=(1,2), plot_title="Moving Minimum Parameters vs Window Size (Gumbell)")


savefig(min_params_1,"Output_Images/parameters_vs_window_size/Moving_min_frechet.pdf")
savefig(min_params_2,"Output_Images/parameters_vs_window_size/Moving_min_gumbell.pdf")
parmparm_block_2[3][1].*log.(block_length)
param_block_av_1 = gumbel_params_av(orbits, window_sizes, x0)
xis_av_1 = pl.plot(block_length, param_block_av_1[1], title=L"ξ", legend=false)
mus_av_1 = pl.plot(block_length, param_block_av_1[2], title=L"μ", legend=false)
theta_av_1 = pl.plot(block_length, param_block_av_1[3], title=L"σ", legend=false)
av_params_1 = pl.plot(xis_av_1, mus_av_1, theta_av_1, size=(1000, 500), layout=(1,3), plot_title="Moving Average Parameters vs Window Size (Frechet)")

param_block_av_2 = params_obs_2_av(orbits, window_sizes)

xis_av_2 = pl.plot(block_length, param_block_av_2[1], title=L"ξ", legend=false)
mus_av_2 = pl.plot(block_length, param_block_av_2[2], title=L"μ", legend=false)
theta_av_2 = pl.plot(block_length, param_block_av_2[3], title=L"σ", legend=false)
av_params_2 = pl.plot(xis_av_2, mus_av_2, theta_av_2, size=(1000, 500), layout=(1,3),  plot_title="Moving Average Parameters vs Window Size (Gumbell)")


savefig(av_params_1,"Output_Images/parameters_vs_window_size/Moving_Average_frechet.pdf")
savefig(av_params_2,"Output_Images/parameters_vs_window_size/Moving_Average_gumbell.pdf")

