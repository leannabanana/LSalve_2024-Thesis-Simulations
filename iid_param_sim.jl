include("chaotic_system_methods.jl")


Random.seed!(1234)

function exponential_distributions(n_vectors, size, distribution)
    random_vectors = [rand(distribution, vector_size) for _ in 1:num_vectors]
    return random_vectors
end

num_vectors = 10
vector_size = 10
λ = 1.0
distribution = Exponential(λ)
X_n = exponential_distributions(num_vectors, vector_size, distribution)
window_sizes = collect(1:4)
function rv_minimum(random_variables, window_sizes)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    for windows in window_sizes
        minimums = map(moving_minimum, random_variables, windows)
        max_min = maximum.(minimums)
        
        shapes = 0.0
        location_1 = 0.0
        scale_1 = 0.0
        
        try
            fit = gevfit(max_min)
            shapes = shape(fit)
            location_1 = location(fit)
            scale_1 = scale(fit)
        catch e
            if isa(e, DomainError)
                shapes = 0.0
                location_1 = 0.0
                scale_1 = 0.0
            else
                rethrow(e)
            end
        end

        append!(shape_params, shapes)
        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
    return shape_params, location_params, scale_params
end

rv_minimum(X_n, window_sizes)