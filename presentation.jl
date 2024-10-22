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
n_iters = 2      # Number of iterations

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
plot(f(x) distances, marker=:o, lw=2, xlabel="Iteration", ylabel="Distance", yscale=:log10, 
     title="Distance Growth between Perturbed and Fixed Points", legend=false)

# Print the initial and final distances to see the expansion
println("Initial distance: ", epsilon)
println("Final distance after ", n_iters, " iterations: ", distances[end])

# Define the map f(x) = 2x mod 1
function f(x)
    return 2 * x % 1 + Normal(0, 1/10^3)
end

# Parameters
n_points = 2  # Number of different initial conditions
n_iters = 2   # Number of iterations

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




function mod_map(x)
    return 2 * x % 1
end
        
function simulate_trajectory(x0, n)
    trajectory = Float64[]
    x = x0
    for i in 1:n
        push!(trajectory, x)
        x = mod_map(x)
    end
    return trajectory
end

# Parameters
n = 1000  # Number of iterations
x0_1 = 0.1  # Initial point for the first trajectory
x0_2 = 0.12  # Initial point for the second trajectory

# Simulate the trajectories
trajectory_1 = simulate_trajectory(x0_1, n)
trajectory_2 = simulate_trajectory(x0_2, n)
trajectory_1
# Plot the trajectories
new = scatter(trajectory_1[1:end-1], trajectory_1[2:end], label = "x0 = $x0_1", mc="red", xlims=(0,1), ylims=(0,1))
scatter!(trajectory_2[1:end-1], trajectory_1[2:end], label = "x0 = $x0_2", mc="pink")


anim = @animate for i in 2:n
    scatter(trajectory_1[1:i-1], trajectory_1[2:i], label = "x0 = $x0_1", mc="red", xlims=(0,1), ylims=(0,1))
    scatter!(trajectory_2[1:i-1], trajectory_2[2:i], label = "x0 = $x0_2", mc="pink")
    xlabel!(L"x_n")
    ylabel!(L"x_{n+1}")
    title!(L"T(x) = 2x(\text{mod}1)")
end

gif(anim, "chaos_animation.gif", fps = 2)

xlabel!(L"x_n")
ylabel!(L"x_{n+1}")
title!(L"T(x) = 2x(mod1)")

savefig(new, "new_chaos.pdf")


# Define the parameters for the GEV distributions
mu = 0          # location parameter
sigma = 1       # scale parameter

# Create the three types of GEV distributions
gumbel = GeneralizedExtremeValue(0, 1, 0)          # Gumbel (Type I)
frechet = GeneralizedExtremeValue(0, 1, 1/2)        # Fréchet (Type II)
weibull = GeneralizedExtremeValue(0, 1, -1/2)        # Weibull (Type III)

# Generate x values
x = -5:0.1:5  # Range for x values

# Calculate PDFs
pdf_gumbel = pdf(gumbel, x)
pdf_frechet = pdf(frechet, x)
pdf_weibull = pdf(weibull, x)

# Plot the PDFs
gev = plot(x, pdf_gumbel, label="Gumbel; ξ = 0", color=:royalblue)
plot!(x, pdf_frechet, label="Fréchet;  ξ = 1", color=:hotpink)
plot!(x, pdf_weibull, label="Weibull;  ξ = -1", color=:purple)

# Add titles and labels
title!("PDFs of GEV Distributions")
xlabel!(L"x")
ylabel!(L"f(x)")


savefig(gev,"Output_Images/18-08-2024/gumbel_min_1sqrt2.pdf")

# Show the plot
display(plot)