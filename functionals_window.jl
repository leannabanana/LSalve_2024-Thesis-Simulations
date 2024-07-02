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
window_size = 14


x0s = rand()

### Simulate 1000 orbits
orbits = simulate_orbits(n_orbits, a, initial_value, interations, pertubation)

###### Frechet Distribution
observable_values = map(orbit -> observable_one(orbit, x0, alpha), orbits) ## Comopse all orbits/trajectories by our observable

#moving minimum functional
frechet_moving_min = maximum.(moving_minimum.(observable_values, window_size))
gev_max_av = maximum_values(frechet_moving_min, 50)
frechet_av = gevfit(gev_max_av)

#moving average functional
frechet_moving_av = maximum.(moving_average.(observable_values, window_size))
gev_max_av = maximum_values(frechet_moving_av, 50)
frechet_av = gevfit(gev_max_av)

###### Gumbell Distribution
observable_values_2 = map(orbit -> observable_two(orbit, x0), orbits)

#moving minimum functional
mov_min_2 = maximum.(moving_minimum.(observable_values_2, 100))
gev_max_2 = maximum_values(mov_min_2, 50)
gumbel_min = gevfit(gev_max_2)

#moving average functional
max_min_av_2 = maximum.(moving_average.(observable_values_2, window_size))
gev_max_av_2 = maximum_values(max_min_av_2, 50)
gumbel_av = gevfit(gev_max_av_2)


#d4 = diagnosticplots(frechet_moving_min)
#d5 = diagnosticplots(frechet_moving_av)
#d6 = diagnosticplots(gumbel_min)
#d7 = diagnosticplots(gumbel_av)

#Saving Diagnostic tests
#draw(PDF("Output_Images/gev_diagnostic_tests/frechet_moving_min"*Date*".pdf", 25cm, 15cm), d4)
#draw(PDF("Output_Images/gev_diagnostic_tests/frechet_moving_av"*Date*".pdf", 25cm, 15cm), d5)
#draw(PDF("Output_Images/gev_diagnostic_tests/gumbell_moving_min"*Date*".pdf",25cm, 15cm), d6)
#draw(PDF("Output_Images/gev_diagnostic_tests/gumbell_moving_av"*Date*".pdf", 25cm, 15cm), d7)

### Observable 3 CURRENTLY BROKEN!!!!
#observable_values_3 = map(orbit -> observable_three(orbit, x0, a, alpha), orbits)
#max_min_3 = maximum.(moving_minimum.(observable_values_3, 90))
#gev_max_3 = maximum_values(max_min_3, 50)
#gevfit(gev_max_3)

#max_min_av_3 = maximum.(moving_average.(observable_values_3, 60))
#gev_max_av_3 = maximum_values(max_min_av_3, 50)
#gevfit(gev_max_av_3)
