
include("chaotic_system_methods.jl")
Random.seed!(1234)

### Define the window sizes
window_sizes = collect(1:50)

# Variables
orbit_length = 10^4
n_orbits = 10^3
block_length = orbit_length ./ window_sizes
initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)

# Parameters
a = 2 
perturbation = 1/10^3

# Simulate orbits
orbits = simulate_orbits(initial_conditions, a, orbit_length, perturbation, n_orbits)

##### Moving Minimums
#Observable 1
function params_min_1(orbits, window_sizes)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    for windows in window_sizes
        observable_values = map(orbit -> observable_one(orbit, 0, 1/3), orbits) ## Compose all orbits/trajectories by our observable
        mov_min_1 = moving_minimum.(observable_values, windows)
        max_min = maximum.(mov_min_1)
        
        shapes = 0.0
        location_1 = 0.0
        scale_1 = 0.0
        
        try
            fit = gevfit(max_min)
            shapes = shape(fit)
            location_1 = location(fit)
            scale_1 = scale(fit)
        catch e
            if isa(e, DomainError)
                shapes = 0.0
                location_1 = 0.0
                scale_1 = 0.0
            else
                rethrow(e)
            end
        end

        append!(shape_params, shapes)
        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
    return shape_params, location_params, scale_params
end

param_block_1 = params_min_1(orbits, window_sizes)

xis = pl.plot(block_length, param_block_1[1], title=L"ξ", legend=false)
mus = pl.plot(block_length, param_block_1[2], title=L"μ", legend=false)
theta = pl.plot(block_length, param_block_1[3], title=L"σ", legend=false)
min_params_1 = pl.plot(xis, mus, theta, size=(1000, 500), layout=(1,3), plot_title="Moving Minimum Parameters vs Window Size (Frechet)")


# Observable 2
function params_obs_2(orbits, window_sizes)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    for windows in window_sizes
        observable_values = map(orbit -> observable_two(orbit, 0), orbits) ## Compose all orbits/trajectories by our observable
        mov_min_1 = moving_minimum.(observable_values, windows)
        gev_max = maximum.(mov_min_1)

        shapes = 0.0
        location_1 = 0.0
        scale_1 = 0.0
        
        try
            fit = gumbelfit(gev_max)
            shapes = shape(fit)
            location_1 = location(fit)
            scale_1 = scale(fit)
        catch e
            if isa(e, DomainError)
                shapes = 0.0
                location_1 = 0.0
                scale_1 = 0.0
            else
                rethrow(e)
            end
        end

        append!(shape_params, shapes)
        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
    return shape_params, location_params, scale_params
end

parm_block_2 = params_obs_2(orbits, window_sizes)
xis_2 = pl.plot(block_length, parm_block_2[1], title=L"ξ", legend=false)
mus_2 = pl.plot(block_length, parm_block_2[2], title=L"μ", legend=false)
theta_2 = pl.plot(block_length, parm_block_2[3], title=L"θ", legend=false)
min_params_2 = pl.plot(xis_2, mus_2, theta_2, size=(1000, 500), layout=(1,3), plot_title="Moving Minimum Parameters vs Window Size (Gumbell)")

savefig(min_params_1,"Output_Images/parameters_vs_window_size/Moving_min_frechet.pdf")
savefig(min_params_2,"Output_Images/parameters_vs_window_size/Moving_min_gumbell.pdf")

### Moving Average
# Observable 1
function params_obs_1_av(orbits, window_size)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    for windows in window_sizes
        observable_values = map(orbit -> observable_one(orbit, 0, α), orbits) ## Comopse all orbits/trajectories by our observable
        mov_min_1= (moving_average.(observable_values, windows))
        max_min = maximum.(mov_min_1)

        shapes = shape(gevfit(max_min))
        location_1 = location(gevfit(max_min))
        scale_1 = scale(gevfit(max_min))
        append!(shape_params, shapes)
        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
    return shape_params, location_params, scale_params
end

param_block_av_1 = params_obs_1_av(orbits, window_sizes)
xis_av_1 = pl.plot(block_length, param_block_av_1[1], title=L"ξ", legend=false)
mus_av_1 = pl.plot(block_length, param_block_av_1[2], title=L"μ", legend=false)
theta_av_1 = pl.plot(block_length, param_block_av_1[3], title=L"σ", legend=false)
av_params_1 = pl.plot(xis_av_1, mus_av_1, theta_av_1, size=(1000, 500), layout=(1,3), plot_title="Moving Average Parameters vs Window Size (Frechet)")

#Observable 2
function params_obs_2_av(orbits, window_size)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]

    for windows in window_sizes
        observable_values = map(orbit -> observable_two(orbit, 0), orbits) ## Comopse all orbits/trajectories by our observable
        mov_min_1= (moving_average.(observable_values, windows))
        max_min = maximum.(mov_min_1)

        shapes = shape(gumbelfit(max_min))
        location_1 = location(gevfit(max_min))
        scale_1 = scale(gevfit(max_min))

        append!(shape_params, shapes)
        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
    return shape_params, location_params, scale_params
end

param_block_av_2 = params_obs_2_av(orbits, window_sizes)

xis_av_2 = pl.plot(block_length, param_block_av_2[1], title=L"ξ", legend=false)
mus_av_2 = pl.plot(block_length, param_block_av_2[2], title=L"μ", legend=false)
theta_av_2 = pl.plot(block_length, param_block_av_2[3], title=L"σ", legend=false)
av_params_2 = pl.plot(xis_av_2, mus_av_2, theta_av_2, size=(1000, 500), layout=(1,3),  plot_title="Moving Average Parameters vs Window Size (Gumbell)")


savefig(av_params_1,"Output_Images/parameters_vs_window_size/Moving_Average_frechet.pdf")
savefig(av_params_2,"Output_Images/parameters_vs_window_size/Moving_Average_gumbell.pdf")

