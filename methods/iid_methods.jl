include("chaotic_system_methods.jl")


#Simulates Random variables
function generate_rv(n_vectors, vector_size, distribution)
    random_vectors = [rand(distribution, vector_size) for _ in 1:n_vectors]
    return random_vectors
end

# IId case applies a moving average and or moving minumum. This uses a MATRIX

function RV_average(random_variables, window_size)
    # location_params = Vector{Float64}(undef,19)
    # scale_params = Vector{Float64}(undef, 19)

    location_params = Float64[]
    scale_params = Float64[]

    for i in 10:size(random_variables)[1]
        minimums = moving_average_matrix(random_variables[1:i, :], window_size)
        max_min = maximum(minimums, dims=1)[:]

        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
     return location_params, scale_params
end

function RV_minimum(random_variables, window_size)
    # location_params = Vector{Float64}(undef,19)
    # scale_params = Vector{Float64}(undef, 19)

    location_params = Float64[]
    scale_params = Float64[]

    for i in 10:size(random_variables)[1]
        minimums = moving_minimum_matrix(random_variables[1:i, :], window_size)
        max_min = maximum(minimums, dims=1)[:]

        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
     return location_params, scale_params
end

##### LIST OF LISTS

function list_average(random_variables, window_size)
    # location_params = Vector{Float64}(undef, 100)
    # scale_params = Vector{Float64}(undef, 100)

    location_params = Float64[]
    scale_params = Float64[]

    Threads.@threads for i in 2:length(random_variables)
        minimums = moving_minimum.(random_variables[1:i], window_size)
        max_min = maximum.(minimums)

        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)
    
    end
    return location_params, scale_params
end

function list_min(random_variables, window_size)
    # location_params = Vector{Float64}(undef,19)
    # scale_params = Vector{Float64}(undef, 19)

    location_params = Float64[]
    scale_params = Float64[]

    Threads.@threads for i in 10:length(random_variables[1])
        data = [variable[1:i] for variable in random_variables]
        minimums = @. moving_minimum(data, window_size)
        max_min = @. maximum(minimums)

        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)

    end
     return location_params, scale_params
end

