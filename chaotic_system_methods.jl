"""
This file initialises methods required to simulate chaotic systems
"""

## Plots and Gadly are both plotting packages a plot command

import Plots as pl #The Extremes package uses Gadfly as its plotting package this distinguishes packaes
using Plots, Extremes, Distributions, Gadfly, LaTeXStrings, Cairo, Fontconfig, Random

set_default_plot_size(25cm, 20cm) ### Choosing a default plot size

##### Simulating our Chaotic Map 
function chaotic_map(a, n_steps) #gives normal map
    y_values = Float64[]
    
    for i in 0:1/n_steps:1
        x = (a*i) % 1
        push!(y_values, x)

    end
    return y_values
end

function T_mod_map(a, start, n) #gives composition n times
    y_values = Float64[]
    x = start

    for i in 0:1/n:1
        x = (a*x) % 1
        push!(y_values, x)

    end
    return y_values
end

### Perturbed System
function T_mod_map_noisey(a, n, pertubation)  #gives composition n times with perturbations
    y_values = Float64[]
    normal_dist_t = Normal(0, pertubation)
    x = 1/n

    for i in 1:n
        x = (a*x) % 1 + rand(normal_dist_t)
        push!(y_values, x)

    end
    return y_values
end

#Defining observables
function observable_one(x_values, x0, alpha)
    distances = abs.(x_values .- x0).^(-alpha)
    return distances
end

function observable_two(x_values, x0)
    distances = -log.(abs.(x_values .- x0))
    return distances
end

function observable_three(x_values, x0, a, alpha)
    distances = a.- (abs.(x_values .- x0).^(-alpha))
    return distances
end