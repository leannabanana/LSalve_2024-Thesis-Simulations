
include("chaotic_system_methods.jl")
include("GEV_fitting.jl")
include("functionals_window.jl")

### Define the window sizes
window_sizes = collect(13:50)

### Simulate 1000 orbits
orbits = simulate_orbits(n_orbits, a, initial_value, interations, pertubation)
#### For whatever reason scale is not defined 
####Observable 1

function params_obs_1(orbits, window_size)
    shape_params = Float64[]
    location_params = Float64[]

    for windows in window_sizes
        observable_values = map(orbit -> observable_one(orbit, x0, alpha), orbits) ## Comopse all orbits/trajectories by our observable
        mov_min_1= (moving_minimum.(observable_values, windows))
        max_min = maximum.(mov_min_1)
        gev_max = maximum_values(max_min, 50)
        shapes = shape(gevfit(gev_max))
        location_1 = location(gevfit(gev_max))
        #scale_1 = scale(gevfit(gev_max))
        append!(shape_params, shapes)
        append!(location_params, location_1)
        #append!(scales_params, scale_1)
    end
    return shape_params, location_params
end
param_block_1 = params_obs_1(orbits, window_sizes)

xis = pl.plot(window_sizes, param_block_1[1], ylabel=L"ξ", xlabel="Window Size", legend=false)
mus = pl.plot(window_sizes, param_block_1[2], ylabel=L"μ", xlabel="Window Size", legend=false)
window_sizes_params = pl.plot(xis, mus, size=(700, 400))

### Observable 2
function params_obs_2(orbits, window_size)
    shape_params = Float64[]
    location_params = Float64[]
    #scales_params = Float64[]
    for windows in window_sizes
        observable_values = map(orbit -> observable_two(orbit, x0), orbits) ## Comopse all orbits/trajectories by our observable
        mov_min_1= (moving_minimum.(observable_values, windows))
        max_min = maximum.(mov_min_1)
        gev_max = maximum_values(max_min, 50)
        shapes = shape(gevfit(gev_max))
        location_1 = location(gevfit(gev_max))
        #scale_1 = scale(gevfit(gev_max))
        append!(shape_params, shapes)
        append!(location_params, location_1)
        #append!(scales_params, scale_1)
    end
    return shape_params, location_params
end


parm_block_2 = params_obs_2(orbits, window_sizes)

xis_2 = pl.plot(window_sizes, parm_block_2[1], ylabel=L"ξ", xlabel="Window Size", legend=false)
mus_2 = pl.plot(window_sizes, parm_block_2[2], ylabel=L"μ", xlabel="Window Size", legend=false)
window_sizes_params = pl.plot(xis_2, mus_2, size=(700, 400), plot_title="Moving minimum parameters vs Window Size")


function params_obs_1_av(orbits, window_size)
    shape_params = Float64[]
    location_params = Float64[]

    for windows in window_sizes
        observable_values = map(orbit -> observable_one(orbit, x0, alpha), orbits) ## Comopse all orbits/trajectories by our observable
        mov_min_1= (moving_average.(observable_values, windows))
        max_min = maximum.(mov_min_1)
        gev_max = maximum_values(max_min, 50)
        shapes = shape(gevfit(gev_max))
        location_1 = location(gevfit(gev_max))
        #scale_1 = scale(gevfit(gev_max))
        append!(shape_params, shapes)
        append!(location_params, location_1)
        #append!(scales_params, scale_1)
    end
    return shape_params, location_params
end

param_block_av_1 = params_obs_1_av(orbits, window_sizes)

xis_av_1 = pl.plot(window_sizes, param_block_av_1[1], ylabel=L"ξ", xlabel="Window Size", legend=false)
mus_av_1 = pl.plot(window_sizes, param_block_av_1[2], ylabel=L"μ", xlabel="Window Size", legend=false)
window_sizes_params = pl.plot(xis, mus, size=(700, 400))


function params_obs_1_av(orbits, window_size)
    shape_params = Float64[]
    location_params = Float64[]

    for windows in window_sizes
        observable_values = map(orbit -> observable_one(orbit, x0, alpha), orbits) ## Comopse all orbits/trajectories by our observable
        mov_min_1= (moving_average.(observable_values, windows))
        max_min = maximum.(mov_min_1)
        gev_max = maximum_values(max_min, 50)
        shapes = shape(gevfit(gev_max))
        location_1 = location(gevfit(gev_max))
        #scale_1 = scale(gevfit(gev_max))
        append!(shape_params, shapes)
        append!(location_params, location_1)
        #append!(scales_params, scale_1)
    end
    return shape_params, location_params
end