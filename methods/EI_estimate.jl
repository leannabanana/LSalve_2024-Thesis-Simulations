"""
This file contains the method of the EI Estimate
"""

include("chaotic_system_methods.jl")

Random.seed!(1234)
window_sizes = collect(1:15)

# ### Define initial_conditions of length n_orbits
# n_orbits = 10^3
# orbit_length = 10^4
# initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)

# ### Simulate orbits
# a = 2 
# pertubation = 1/10^3
# p0 = 0
# p1 = 1/(sqrt(2*pi))

# orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
# observables = map(orbit -> observable_two(orbit, p1), orbits)
# mov_min_1 = moving_minimum.(observables, 50)
# gev_max = maximum.(mov_min_1)


### Estimate θ

function extremal_FerroSegers(Y, p)
    # This function computes the extremal index theta by using the
    # method proposed by Ferro-Segers (Ferro, C. A. T., and
    # J. Segers (2003), Inference for clusters of extremes,
    # J. R. Stat. Soc., Ser. B, 65, 545-556.).

    # INPUTS:
    # - Y: a vector containing a univariate time series
    # - p: a quantile value
    # OUTPUTS:
    # - theta: the estimate of the extremal index.
    
    # Extract the threshold u corresponding to the quantile p
    u = quantile(Y, p)
    
    # Compute the exceedances
    Si = findall(y -> y > u, Y)
    
    # Compute the cluster lengths
    Ti = diff(Si)
    
    # Compute the total number of clusters
    N = length(Ti)
    
    # Use the Ferro-Segers formula to extract theta
    theta = 2 * (sum(Ti .- 1))^2 / (N * sum((Ti .- 1) .* (Ti .- 2))) 
    
    return theta
end

# extremal_FerroSegers(gev_max, 0.95)

# # Variables
# orbit_length = 10^4
# n_orbits = 10^3
# initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)
# x0 = 0
# α = 1/3

# # Parameters
# a = 2 
# perturbation = 1/10^3

# ### Simulate 1000 orbits
# orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
# observed_orbits = observable_two.(orbits, x0)
# orbit_matrix = reduce(hcat, observed_orbits)
# size(orbit_matrix)
function EI_estimation_min(orbits, window_size)
    # location_params = Vector{Float64}(undef,19)
    # scale_params = Vector{Float64}(undef, 19)

    location_params = Float64[]
    scale_params = Float64[]
    EI = Float64[]

    for i in 10:size(orbits)[1]
        minimums = moving_minimum_matrix(orbits[1:i, :], window_size)
        max_min = maximum(minimums, dims=1)[:]

        EI_estimate = extremal_FerroSegers(max_min, 0.95)
        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)
        append!(EI, EI_estimate)
    end
     return location_params, scale_params, EI
end

# params = EI_estimation_min(orbit_matrix, 10)


# x_axis_1 =  collect(10:size(orbit_matrix)[1])
# σ_min_orbit =scatter(x_axis_1, params[2], legend=false, xlabel = "Block Length", ylabel=L"\sigma^*", ms=1/2, ma =1/2, mc="indianred1", markerstrokecolor="indianred1", title =L"GEV Maxima sampled from $φ(T^n(x))$ with $x_0 = 0$",  gridcolor=:gray19, gridalpha=1/2)
# μ_min_orbit = scatter(x_axis_1, params[1], legend=false, xlabel = "Block Length", ylabel=L"\mu^*", ms=1/2, ma =1/2, mc="indianred1",  markerstrokecolor="indianred1",  title =L"GEV Maxima sampled from $φ(T^n(x))$ with $x_0 = 0$", gridcolor=:gray19, gridalpha=1/2)
# θ_min = scatter(x_axis_1, params[3], legend=false, xlabel = "Block Length", ylabel=L"θ", ms=1.5, ma =1/3, mc="indianred1",  markerstrokecolor="black",  title =L"GEV Maxima sampled from $φ(T^n(x))$ with $x_0 = 0$", gridcolor=:gray19, gridalpha=1/2)

# savefig(σ_min_orbit,"Output_Images/EI_estimation_orbit/σmin_orbit_EI.pdf")
# savefig(μ_min_orbit,"Output_Images/EI_estimation_orbit/μmin_orbit_EI.pdf")
# savefig(θ_min,"Output_Images/EI_estimation_orbit/EI_estimate_min.pdf")


# orbits_1 = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
# observed_orbits = observable_two.(orbits_1, 1/sqrt(2))
# orbit_matrix_1 = reduce(hcat, observed_orbits)

# params_2 = EI_estimation_min(orbit_matrix_1, 10)
# σ_min1_orbit =scatter(x_axis_1, params_2[2], legend=false, xlabel = "Block Length", ylabel=L"\sigma^*", ms=1/2, ma =1/2, mc="indianred1", markerstrokecolor="indianred1", title =L"GEV Maxima sampled from $φ(T^n(x))$ with $x_0 = \dfrac{1}{\sqrt{2}}$",  gridcolor=:gray19, gridalpha=1/2)
# μ_min1_orbit = scatter(x_axis_1, params_2[1], legend=false, xlabel = "Block Length", ylabel=L"\mu^*", ms=1/2, ma =1/2, mc="indianred1",  markerstrokecolor="indianred1",  title =L"GEV Maxima sampled from $φ(T^n(x))$ with $x_0 = \dfrac{1}{\sqrt{2}}$", gridcolor=:gray19, gridalpha=1/2)
# θ1_min = scatter(x_axis_1, params_2[3], legend=false, xlabel = "Block Length", ylabel=L"θ", ms=1.5, ma =1/3, mc="indianred1",  markerstrokecolor="black",  title =L"GEV Maxima sampled from $φ(T^n(x))$ with $x_0 = \dfrac{1}{\sqrt{2}}$", gridcolor=:gray19, gridalpha=1/2)


# savefig(σ_min1_orbit,"Output_Images/EI_estimation_orbit/σmin_orbit_irrational.pdf")
# savefig(μ_min1_orbit,"Output_Images/EI_estimation_orbit/μmin_orbit_irrational.pdf")
# savefig(θ1_min,"Output_Images/EI_estimation_orbit/EI_estimate_min_irrational.pdf")


function EI_estimation_average(orbits, window_size)
    # location_params = Vector{Float64}(undef,19)
    # scale_params = Vector{Float64}(undef, 19)

    location_params = Float64[]
    scale_params = Float64[]
    EI = Float64[]

    for i in 10:size(orbits)[1]
        minimums = moving_average_matrix(orbits[1:i, :], window_size)
        max_min = maximum(minimums, dims=1)[:]

        EI_estimate = extremal_FerroSegers(max_min, 0.95)
        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)
        append!(EI, EI_estimate)
    end
     return location_params, scale_params, EI
end


# av_orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
# observed_orbits_av = observable_two.(av_orbits, 1/sqrt(2))
# orbit_matrix_2 = reduce(hcat, observed_orbits_av)

# av_par = EI_estimation_average(orbit_matrix_2, 10)

# av_orbits_2 = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
# observed_orbits_av2 = observable_two.(av_orbits_2, 0)
# orbit_matrix_3 = reduce(hcat, observed_orbits_av2)

# av_par_2 = EI_estimation_average(orbit_matrix_3, 10)
