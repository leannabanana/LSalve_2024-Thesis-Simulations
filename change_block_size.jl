include("methods/chaotic_system_methods.jl")
include("methods/parameter_methods.jl")

### Define a function which gets me parameters for changing window sizes for the moving minimum functional 

Random.seed!(1234)
### Simulate a bunch of RV's 

# Function to generate exponential distributions cumulatively
function exponential_distributions_cumulative(num_vectors_list, vector_size, distribution)
    cumulative_vectors = []
    all_vectors = []

    for num_vectors in num_vectors_list
        new_vectors = [rand(distribution, vector_size) for _ in 1:num_vectors]
        append!(cumulative_vectors, new_vectors)
        push!(all_vectors, deepcopy(cumulative_vectors))
    end

    return all_vectors
end

# Function to fit Gumbel distribution to moving minimums
function rv_minimum(random_variables, window_size)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    for vectors in random_variables
        minimums = [moving_minimum(vec, window_size) for vec in vectors]
        max_min = maximum.(minimums)

        fit = gumbelfit(max_min)
        shapes = shape(fit)
        location_1 = location(fit)
        scale_1 = scale(fit)

        push!(shape_params, shapes)
        push!(location_params, location_1)
        push!(scale_params, scale_1)
    end
    return shape_params, location_params, scale_params
end
