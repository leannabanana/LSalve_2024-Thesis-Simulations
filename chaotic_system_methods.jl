"""
This file initialises methods required to simulate chaotic systems
"""
## Plots and Gadly are both plotting packages a plot command
import Plots as pl #The Extremes package uses Gadfly as its plotting package this distinguishes packaes
using Plots, Extremes, Distributions, Gadfly, LaTeXStrings, Fontconfig, Random, DataFrames, CSV, Statistics

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
    y_values = [start]
    normal_dist_t = Normal(0, pertubation)
    x = start

    for i in 1:n-1
        x = (a*x) % 1 + rand(normal_dist_t)
        push!(y_values, x) 

    end
    return y_values
end

#Defining observables
function observable_one(x_values, x0, α)
    distances = abs.(x_values .- x0).^(-α)
    return distances
end

function observable_two(x_values, x0)
    distances = -log.(abs.(x_values .- x0))
    return distances
end

function observable_three(x_values, x0, a, α)
    distances = a.- (abs.(x_values .- x0).^(-α))
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

### GEV fitting, this takes the maximum values given k blocks
function maximum_values(data, k)
    blocks = k_blocks(data, k)
    maxima = [maximum(block) for block in blocks]
    return maxima
end

### This function returns the shape parameter depending on how many blocks we have
function xi_params(data, block_number)
    shape_params = Float64[]
    for n in block_number
        shapes = shape(gevfit(maximum_values(data, n)))
        append!(shape_params, shapes)
    end 

    return shape_params
end

moving_average(data, window_size) = [sum(@view data[i:(i+window_size-1)])/window_size for i in 1:(length(data)-(window_size-1))]

moving_minimum(data, window_size) = [minimum(@view data[i:(i+window_size-1)]) for i in 1:(length(data)-(window_size-1))]

function simulate_orbits(initial_conditions::Vector{Float64}, a, n_length, perturbation, num_orbits)
    all_orbits = Vector{Vector{Float64}}(undef, num_orbits)
    for i in 1:num_orbits
        all_orbits[i] = T_mod_map_noisey(a, initial_conditions[i], n_length, perturbation)
    end
    return all_orbits
end
