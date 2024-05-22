
include("chaotic_system_methods.jl")
include("GEV_fitting.jl")
include("functionals_window.jl")


# Define the window sizes
window_sizes = collect(50:40)

####Observable 1
shapeparams_min_1 = Float64[]
shapeparams_av_1 = Float64[]

function shape_parameters(n_orbits, a, initial_value, interations, pertubation, window_size)
    shapeparams_min_1 = Float64[]
    for windows in  window_sizes
        orbits = simulate_orbits(n_orbits, a, initial_value, interations, pertubation)
        observable_values = map(orbit -> observable_one(orbit, x0, alpha), orbits)
        max_min = maximum.(moving_minimum.(observable_values, windows))
        gev_max = maximum_values(max_min, 50)
        gev_xi = shape(gevfit(gev_max))
        append!(shapeparams_min_1, gev_xi)
    end
    return shapeparams_min_1
end


shape_parameters(n_orbits, a, initial_value, interations, pertubation, window_sizes)