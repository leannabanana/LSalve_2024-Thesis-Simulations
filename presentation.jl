include("methods/EI_estimate.jl")
include("methods/chaotic_system_methods.jl")

# Define the map f(x) = 2x mod 1
function f(x)
    return 2 * x % 1
end

# Initialize parameters
x0 = 1/3          # Initial point near the fixed point
epsilon = 1e-5    # Small perturbation (nearby point)
x1 = x0 + epsilon # Perturbed point
n_iters = 10      # Number of iterations

# Arrays to store iterates and distances
x0_vals = [x0]
x1_vals = [x1]
distances = [abs(x1 - x0)]  # Initial distance

# Perform iterations
for i in 1:n_iters
    x0 = f(x0)
    x1 = f(x1)
    push!(x0_vals, x0)
    push!(x1_vals, x1)
    push!(distances, abs(x1 - x0))
end

# Plot the distance growth over iterations
plot(0:n_iters, distances, marker=:o, lw=2, xlabel="Iteration", ylabel="Distance", yscale=:log10, 
     title="Distance Growth between Perturbed and Fixed Points", legend=false)

# Print the initial and final distances to see the expansion
println("Initial distance: ", epsilon)
println("Final distance after ", n_iters, " iterations: ", distances[end])




using Plots

# Define the map f(x) = 2x mod 1
function f(x)
    return 2 * x % 1 + Normal(0, 1/10^3)
end

# Parameters
n_points = 100  # Number of different initial conditions
n_iters = 100   # Number of iterations

# Arrays to store initial conditions and final values after n_iters
initial_conditions = rand(n_points)  # Randomly generate initial conditions in [0, 1]
final_values = zeros(n_points)

# Iterate over each initial condition and compute its trajectory
for i in 1:n_points
    x = initial_conditions[i]
    for n in 1:n_iters
        x = f(x)
    end
    final_values[i] = x  # Store the final value after n_iters
end
final_values
# Plot final values against initial conditions
scatter(initial_conditions, final_values, xlabel="Initial Condition (x_0)", ylabel="Value after $n_iters (x_n)",
        title="Final Values vs Initial Conditions for T(x) = 2x mod 1", legend=false, label="")
