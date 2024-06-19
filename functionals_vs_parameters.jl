
include("chaotic_system_methods.jl")
include("GEV_fitting.jl")
include("functionals_window.jl")

### Define the window sizes
window_sizes = collect(1:25)

### Simulate 1000 orbits
orbits = simulate_orbits(n_orbits, a, initial_value, interations, pertubation)
#### For whatever reason scale is not defined 


##### Moving Minimums
#Observable 1

function params_obs_1(orbits, window_sizes)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    for windows in window_sizes
        observable_values = map(orbit -> observable_one(orbit, x0, alpha), orbits) ## Compose all orbits/trajectories by our observable
        mov_min_1 = moving_minimum.(observable_values, windows)
        max_min = maximum.(mov_min_1)
        gev_max = maximum_values(max_min, 50)
        
        shapes = 0.0
        location_1 = 0.0
        scale_1 = 0.0
        
        try
            fit = gevfit(gev_max)
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

param_block_1 = params_obs_1(orbits, window_sizes)
param_block_1[3]
xis = pl.plot(window_sizes, param_block_1[1], title=L"ξ", legend=false)
mus = pl.plot(window_sizes, param_block_1[2], title=L"μ", legend=false)
theta = pl.plot(window_sizes, param_block_1[3], title=L"θ", legend=false)
min_params_1 = pl.plot(xis, mus, theta, size=(1000, 500), layout=(1,3), plot_title="Moving Average Parameters vs Window Size (Observable 1)")


# Observable 2

function params_obs_2(orbits, window_sizes)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    for windows in window_sizes
        observable_values = map(orbit -> observable_two(orbit, x0), orbits) ## Compose all orbits/trajectories by our observable
        mov_min_1 = moving_minimum.(observable_values, windows)
        max_min = maximum.(mov_min_1)
        gev_max = maximum_values(max_min, 50)
        
        shapes = 0.0
        location_1 = 0.0
        scale_1 = 0.0
        
        try
            fit = gevfit(gev_max)
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

xis_2 = pl.plot(window_sizes, parm_block_2[1], title=L"ξ", legend=false)
mus_2 = pl.plot(window_sizes, parm_block_2[2], title=L"μ", legend=false)
theta_2 = pl.plot(window_sizes, param_block_1[3], title=L"θ", legend=false)
min_params_2 = pl.plot(xis_2, mus_2, theta_2, size=(1000, 500), layout=(1,3), plot_title="Moving Minimum Parameters vs Window Size (Observable 2)")

savefig(min_params_1,"Output_Images/parameters_vs_window_size/Moving_min_obs1.png")
savefig(min_params_2,"Output_Images/parameters_vs_window_size/Moving_min_obs2.png")

### Moving Average
# Observable 1
function params_obs_1_av(orbits, window_size)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    for windows in window_sizes
        observable_values = map(orbit -> observable_one(orbit, x0, alpha), orbits) ## Comopse all orbits/trajectories by our observable
        mov_min_1= (moving_average.(observable_values, windows))
        max_min = maximum.(mov_min_1)
        gev_max = maximum_values(max_min, 50)
        shapes = shape(gevfit(gev_max))
        location_1 = location(gevfit(gev_max))
        scale_1 = scale(gevfit(gev_max))
        append!(shape_params, shapes)
        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
    return shape_params, location_params, scale_params
end

param_block_av_1 = params_obs_1_av(orbits, window_sizes)
xis_av_1 = pl.plot(window_sizes, param_block_av_1[1], title=L"ξ", legend=false)
mus_av_1 = pl.plot(window_sizes, param_block_av_1[2], title=L"μ", legend=false)
theta_av_1 = pl.plot(window_sizes, param_block_av_1[3], title=L"θ", legend=false)
av_params_1 = pl.plot(xis_av_1, mus_av_1, theta_av_1, size=(1000, 500), layout=(1,3), plot_title="Moving Average Parameters vs Window Size (Observable 1)")

#Observable 2
function params_obs_2_av(orbits, window_size)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]

    for windows in window_sizes
        observable_values = map(orbit -> observable_two(orbit, x0), orbits) ## Comopse all orbits/trajectories by our observable
        mov_min_1= (moving_average.(observable_values, windows))
        max_min = maximum.(mov_min_1)
        gev_max = maximum_values(max_min, 50)
        shapes = shape(gevfit(gev_max))
        location_1 = location(gevfit(gev_max))
        scale_1 = scale(gevfit(gev_max))

        append!(shape_params, shapes)
        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
    return shape_params, location_params, scale_params
end

param_block_av_2 = params_obs_2_av(orbits, window_sizes)

xis_av_2 = pl.plot(window_sizes, param_block_av_2[1], title=L"ξ", legend=false)
annotate!(0, 0.5, text("Y Axis Title", rotation=90, halign=:center, valign=:center))
mus_av_2 = pl.plot(window_sizes, param_block_av_2[2], title=L"μ", legend=false)

theta_av_2 = pl.plot(window_sizes, param_block_av_2[3], title=L"θ", legend=false)
av_params_2 = pl.plot(xis_av_2, mus_av_2, theta_av_2, size=(1000, 500), layout=(1,3),  plot_title="Moving Average Parameters vs Window Sie (Observable 2)")


savefig(av_params_1,"Output_Images/parameters_vs_window_size/Moving_Average_obs1.png")
savefig(av_params_2,"Output_Images/parameters_vs_window_size/Moving_Average_obs2.png")

