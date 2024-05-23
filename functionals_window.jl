"""
This file fits a GEV to functionals
"""

include("chaotic_system_methods.jl")
include("GEV_fitting.jl")

#Define Constants
a = 3
alpha = 1/6
x0 = 1/3
c = 2
interations = 10^3
pertubation = 1/(10^3)
initial_value = 0.01
n_orbits = 10^3
window_size = 30


### Simulate 1000 orbits
orbits = simulate_orbits(n_orbits, a, initial_value, interations, pertubation)

### Observable 1
observable_values = map(orbit -> observable_one(orbit, x0, alpha), orbits) ## Comopse all orbits/trajectories by our observable
mov_min_1= (moving_minimum.(observable_values, 90))
max_min = maximum.(mov_min_1)
gev_max = maximum_values(max_min, 50)
gevfit(gev_max)

max_min_av_1 = maximum.(moving_average.(observable_values, 90))
gev_max_av = maximum_values(max_min_av_1, 50)
gevfit(gev_max_av)

### Observable 2
observable_values_2 = map(orbit -> observable_two(orbit, x0), orbits)
mov_min_2 = maximum.(moving_minimum.(observable_values_2, 100))
gev_max_2 = maximum_values(mov_min_2, 50)
gevfit(gev_max_2)

max_min_av_2 = maximum.(moving_average.(observable_values_2, 90))
gev_max_av_2 = maximum_values(max_min_av_2, 50)
gevfit(gev_max_av_2)

### Observable 3 CURRENTLY BROKEN!!!!
#observable_values_3 = map(orbit -> observable_three(orbit, x0, a, alpha), orbits)
#max_min_3 = maximum.(moving_minimum.(observable_values_3, 90))
#gev_max_3 = maximum_values(max_min_3, 50)
#gevfit(gev_max_3)

#max_min_av_3 = maximum.(moving_average.(observable_values_3, 60))
#gev_max_av_3 = maximum_values(max_min_av_3, 50)
#gevfit(gev_max_av_3)
