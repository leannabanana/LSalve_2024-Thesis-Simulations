"""
This file contains the method of the EI Estimate
"""

include("chaotic_system_methods.jl")

Random.seed!(1234)

### Define initial_conditions of length n_orbits
n_orbits = 100
orbit_length = 100
initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)

### Simulate orbits
a = 2 
pertubation = 1/10^3
p0 = 0
p1 = 1/(sqrt(2))

orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
observables = map(orbit -> observable_two(orbit, p1), orbits)
mov_min_1 = moving_minimum.(observables, 15)
gev_max = maximum.(mov_min_1)


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

extremal_FerroSegers(gev_max, 0.95)

function EI_estimate_plot(orbits, window_size, x0)
    EI = Float64[]
    location = Float64[]
    scale = Float64[]

    Threads.@threads for i in 2:length(orbits)
        observable = map(orbit -> observable_two(orbit, x0), orbits)
        min = moving_minimum.(observable, window_size)
        max = maximum.(min)

        EI_estimate = extremal_FerroSegers(max, 0.95)
        fit = gumbelfit(max)
        location_pm = location(fit)
        scale_pm = scale(fit)
    
        append!(location, location_pm)
        append!(scale, scale_pm)
        append!(EI, EI_estimate)
    end
    return EI, location, scale
end

# @btime EI_estimate_plot(orbits, 10, 0)
println("x0 = 0 gets us 0.776479812607765 and x0 = 1/sqrt(2) = 1.0181")

function EI_window(orbits, window_size, x0)
    EI = Float64[]
    location = Float64[]
    scale = Float64[]

    Threads.@threads for windows in window_sizes
        observable = map(orbit -> observable_two(orbit, x0), orbits)
        min = moving_minimum.(observable, windows)
        max = maximum.(min)

        EI_estimate = extremal_FerroSegers(max, 0.95)
        fit = gumbelfit(max)
        location_pm = location(fit)
        scale_pm = scale(fit)
    
        append!(location, location_pm)
        append!(scale, scale_pm)
        append!(EI, EI_estimate)
    end
    return EI, location, scale
end


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