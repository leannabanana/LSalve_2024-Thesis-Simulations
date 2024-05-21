"""
This file simulates 1000 orbits of length 10^3
"""

include("chaotic_system_methods.jl")
#Define Constants
a = 3
alpha = 1/6
x0 = 1/3
c = 2
interations = 10^3
pertubation = 1/(10^3)
initial_value = 0.01
n_orbits = 10^3

orbits = simulate_orbits(n_orbits, a, initial_value, interations, pertubation)

moving_minimum(orbits, 400)


gevfit(maximum_values(moving_average(orbits[4], 400), 400))