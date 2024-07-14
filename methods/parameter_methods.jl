include("chaotic_system_methods.jl")

##### Moving Minimums
#Frechet Observable
function frechet_params_min(orbits, window_sizes, x0, α)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    for windows in window_sizes
        observable_values = map(orbit -> observable_one(orbit, x0, α), orbits) ## Compose all orbits/trajectories by our observable
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


# Gumbel Observable
function gumbel_params_min(orbits, window_sizes, x0)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    for windows in window_sizes
        observable_values = map(orbit -> observable_two(orbit, x0), orbits) ## Compose all orbits/trajectories by our observable
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

### Moving Average
#Frechet Observable
function frechet_params_av(orbits, window_sizes, x0, α)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    for windows in window_sizes
        observable_values = map(orbit -> observable_one(orbit, x0, α), orbits) ## Comopse all orbits/trajectories by our observable
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


#Gumbel Observable
function gumbel_params_av(orbits, window_size, x0)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]

    for windows in window_sizes
        observable_values = map(orbit -> observable_two(orbit, x0), orbits) ## Comopse all orbits/trajectories by our observable
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

