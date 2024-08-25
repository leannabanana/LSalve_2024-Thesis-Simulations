import Plots as pl #The Extremes package uses Gadfly as its plotting package this distinguishes packaes
using Plots, Extremes, Distributions, Random, DataFrames, CSV, Statistics, Base.Threads

######## Functions

function T_mod_map_noisey(a, start, n, pertubation)  #gives composition n times with perturbations
    y_values = [start]
    normal_dist_t = Normal(0, pertubation)
    x = start

    for i in 1:n
        pain = rand(normal_dist_t)
        x = (a*x) % 1 + pain
        push!(y_values, x) 
        
    end
    return y_values
end

function observable_two(x_values, x0)
    distances = -log.(abs.(x_values .- x0))
    return distances
end


function extremal_FerroSegers(Y, p)
    # This function computes the extremal index theta by using the
    # method proposed by Ferro-Segers (Ferro, C. A. T., and
    # J. Segers (2003), Inference for clusters of extremes,
    # J. R. Stat. Soc., Ser. B, 65, 545-556.).

    # INPUTS:
    # - Y: a vector containing a univariate time series
    # - p: a quantile value
    # OUTPUTS:
    # - theta: the estimate of the extremal index.
    
    # Extract the threshold u corresponding to the quantile p
    u = quantile(Y, p)
    
    # Compute the exceedances
    Si = findall(y -> y > u, Y)
    
    # Compute the cluster lengths
    Ti = diff(Si)
    
    # Compute the total number of clusters
    N = length(Ti)
    
    # Use the Ferro-Segers formula to extract theta
    theta = 2 * (sum(Ti .- 1))^2 / (N * sum((Ti .- 1) .* (Ti .- 2))) 
    
    return theta
end


function moving_average_matrix(data::Matrix{Float64}, window_size::Int)
    # Get the number of columns and rows
    n_rows, n_cols = size(data)
    
    # Create an empty matrix to store the result
    result = Matrix{Float64}(undef, n_rows - window_size + 1, n_cols)
    
    # Apply moving minimum to each column
    for j in 1:n_cols
        result[:, j] = moving_average(view(data[:, j], :), window_size)
    end
    
    return result
end


function simulate_orbits(initial_conditions::Vector{Float64}, a, n_length, perturbation, num_orbits)
    all_orbits = Vector{Vector{Float64}}(undef, num_orbits)
    for i in 1:num_orbits
        all_orbits[i] = T_mod_map_noisey(a, initial_conditions[i], n_length, perturbation)
    end
    return all_orbits
end

function moving_minimum_matrix(data::Matrix{Float64}, window_size::Int)
    # Get the number of columns and rows
    n_rows, n_cols = size(data)
    
    # Create an empty matrix to store the result
    result = Matrix{Float64}(undef, n_rows - window_size + 1, n_cols)
    
    # Apply moving minimum to each column
    for j in 1:n_cols
        result[:, j] = moving_minimum(view(data[:, j], :), window_size)
    end
    
    return result
end

moving_minimum(data, window_size) = [minimum(@view data[i:(i+window_size-1)]) for i in 1:(length(data)-(window_size-1))]


function f_EI_window_min(orbits, window_sizes)
    location_params = Vector{Vector{Float64}}(undef, length(window_sizes))
    scale_params = Vector{Vector{Float64}}(undef, length(window_sizes))
    EI = Vector{Float64}(undef, length(window_sizes))
    shape_params = Vector{Vector{Float64}}(undef, length(window_sizes))
    
    Threads.@threads for i in 1:length(window_sizes)
        windows = window_sizes[i]
        minimums = moving_minimum_matrix(orbits, windows)
        max_min = maximum(minimums, dims=1)[:]

        fit = gumbelfit(max_min)

        EI_estimate = extremal_FerroSegers(max_min, 0.95)
        location_1 = location(fit)
        scale_1 = scale(fit)
        shapes = shape(fit)

        location_params[i] = location_1
        scale_params[i] = scale_1
        EI[i] = EI_estimate
        shape_params[i] = shapes
    end
    
    return location_params, scale_params, EI, shape_params
end



######## The actual code doing things

Random.seed!(1234)

n_orbits = 10^3
orbit_length = 10^5
initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)

window_sizes = collect(1:10)
a = 2
pertubation = 1/10^3
p0 = 0
p1 = 1/(sqrt(2))

orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
average_orbits = observable_two.(orbits, p1)
mat_orb3 = reduce(hcat, average_orbits)


testing = f_EI_window_min(mat_orb3, window_sizes)


parameters = DataFrame(location = vcat(testing[1]...), scale =  vcat(testing[2]...), EI = vcat(testing[3]...))

CSV.write("gumbel_data_9.csv", parameters)

df0 = CSV.read("gumbel_data_1.csv", DataFrame)
df0.scale

scatter(window_sizes, df0.scale)
pl.plot!(window_sizes, df0.scale[1] ./ 2 .^(window_sizes .- 1))

scatter(window_sizes, df0.location)
pl.plot!(window_sizes, df0.location[1] ./ window_sizes )