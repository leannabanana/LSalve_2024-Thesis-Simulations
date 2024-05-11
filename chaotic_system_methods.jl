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
    start = 1/n_steps
    
    for i in start:1/n_steps:1
        x = (a*i) % 1
        push!(y_values, x)

    end
    return y_values
end

function T_mod_map(a, start, n) #gives composition n times
    y_values = Float64[]
    x = start
    
    for i in 1:n
        x = (a*x) % 1
        push!(y_values, x) 

    end
    return y_values
end

### Perturbed System
function T_mod_map_noisey(a, start, n, pertubation)  #gives composition n times with perturbations
    y_values = Float64[]
    normal_dist_t = Normal(0, pertubation)
    x = start

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

function k_blocks(data, k)
    n = length(data)
    q, r = divrem(n, k)  # Calculate the quotient and remainder
    blocks = []
    start = 1
    for i = 1:k
        len = q + (i <= r ? 1 : 0)  # Blocks indexed <= r have one more element
        push!(blocks, data[start:(start+len-1)])
        start += len
    end
    blocks
end

function maximum_values(data, k)
    blocks = k_blocks(data, k)
    maxima = [maximum(block) for block in blocks]
    return maxima
end