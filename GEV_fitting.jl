"""
This file outputs simulations + fits GEV according to methods in "chaotic_system_methods.jl"
"""

include("methods/chaotic_system_methods.jl")
Date = "_22-06-24_"
Random.seed!(1234)


#Define Constants
a = 2
α = 1/3
x0 = 1/3
c = 2
interations = 10^3
pertubation = 1/(10^3)
initial_value = 1/sqrt(2)

#Create T(x) = 3xmod1
y_map = chaotic_map(a, interations)
x_map = range(0, 1, length=length(y_map))

### T^n(x)
y_n = T_mod_map(a, initial_value, interations)
y_n_noise = T_mod_map_noisey(a, initial_value,  interations, pertubation)

### Perturbed Chaotic Map 3x mod 1 with our 3 observables
y_obs_1 = observable_one(T_mod_map_noisey(a, initial_value, interations, pertubation), x0 , α)
y_obs_2 = observable_two(T_mod_map_noisey(a,  initial_value, interations, pertubation), x0 )
y_obs_3 = observable_three(T_mod_map_noisey(a,  initial_value, interations, pertubation), x0 , c, α)


y_n_noise = T_mod_map_noisey(a, 1/3,  100000, 0)

### Plit into k_blocks and find its maximum_values
maximums_1 = maximum_values(y_obs_1, 50)
maximums_2 = maximum_values(y_obs_2, 50)
maximums_3 = maximum_values(y_obs_3, 50)



test = T_mod_map_noisey(2, 1/π, 10, 0)
### We plot the scatter plots onto T(x)
#Observable 1
orbit_1 = scatter(y_n_noise[1:end-1], y_n_noise[2:end], ms=2, ma=0.5, mc=:lightblue,  xlims=(0,1), legend=false)

p1 = scatter!(x_map, test, mc =:pink, ms=2, ma=0.5, title=L"\phi = d(x, x_0)^{-α}")

#Observable 2
orbit_2 = scatter(y_n_noise[1:end-1] ,y_n_noise[2:end], ms=2, ma=0.5, mc=:lightblue, xlims=(0,1), legend=false)
p2 = scatter!(x_map, y_obs_2, mc =:pink, ms=2, ma=0.5, title=L"\phi = -log(d(x, x_0))")

#Observable 3
orbit_3 = scatter(y_n_noise[1:end-1],y_n_noise[2:end], ms=2, ma=0.5, mc=:lightblue, xlims=(0,1), legend=false)
p3 = scatter!(x_map, y_obs_3, mc =:pink, ms=2, ma=0.5, title=L"\phi = c - d(x, x_0)^{-α}")

### Put 3 plots onto one image 
orbiting = pl.plot(orbit_1, orbit_2, orbit_3, layout=(1,3), size=(1100,500))

#savefig(orbiting,"Output_Images/observable_scatter_plots/perturbed_observables.pdf")

### Fitting a GEV
#Observable 1
fit_obs1 = gevfit(maximums_1)
d11 = histplot(fit_obs1)
d1 = diagnosticplots(fit_obs1)

shape(fit_obs1)
scale(fit_obs1)
#Observable 2
fit_obs2 = gevfit(maximums_2)
d2 = diagnosticplots(fit_obs2)

#Observable 3
fit_obs3 = gevfit(maximums_3)
d3 = diagnosticplots(fit_obs3)

#Save Diagnostic Tests
draw(PDF("Output_Images/gev_diagnostic_tests/frecheee.pdf",25cm, 15cm), d11)
#draw(PDF("Output_Images/gev_diagnostic_tests/gumbell"*Date*".pdf",25cm, 15cm), d2)
#draw(PDF("Output_Images/gev_diagnostic_tests/weibull"*Date*".pdf", 25cm, 15cm), d3)
Date = 3
"Output_Images/gev_diagnostic_tests/gumbell"* string(Date) *".pdf"


function doubling_map(x::Float64, n::Int)
    trajectory = Float64[] # Store the results
    for i in 1:n
        x = 2x % 1
        push!(trajectory, x)
    end
    return trajectory
end

x0 = 1/3 # Initial value
n = 20  # Number of iterations

trajectory = doubling_map(1/3, n)

scatter(trajectory[1:3-1], trajectory[2:3], legend=false, mc= "lightcoral", xlabel = "n (iterations)", ylabel = "2x mod1", xlims=(0,1), ylims=(0,1))
trajectory[1:3-1]
trajectory[2:3]

anim = @animate for i ∈ 1:3
    scatter(trajectory[1:i-1], trajectory[2:i], xlims=(0,1), ylims=(0,1), legend=false, mc= "lightcoral", xlabel = "n (iterations)", ylabel = "2x mod1")
end

gif(anim, "2xmod1_0_1.gif", fps = 2)


# Doubling map: 2x mod 1
function doubling_map(x)
    return 2 * x % 1
end

# Parameters
n_iter = 100  # Number of iterations
x0_1 = 1/3   # Initial condition 1
x0_2 = 0.1001 # Initial condition 2 (slightly different)

# Initialize arrays to store results
x1 = zeros(n_iter)
x2 = zeros(n_iter)
x1[1] = x0_1
x2[1] = x0_2

# Iterate the map
for i in 2:n_iter
    x1[i] = doubling_map(x1[i-1])
    x2[i] = doubling_map(x2[i-1])
end


n_orbits = 2
orbit_length =  2

initial_conditions = rand(n_orbits)
# initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)
a = 2
pertubation = 1/10^3
p0 = 0
p1 = 1/(sqrt(2))
orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
mat_orb3 = reduce(hcat, orbits)
mat_orb3
pog = scatter(mat_orb3[1:end-1, 1], mat_orb3[2:end, 1], label="Trajectory 1", color=:blue, markersize=3, markerstrokewidth =0.5, legend=false, xlim=(0,1), ylim=(0,1), xlabel=L"x_{n}", ylabel=L"x_{n+1}")

for i in 2:n_orbits
    scatter!(mat_orb3[1:end-1, i], mat_orb3[2:end, i], label="Trajectory $i", markersize=3, markerstrokewidth =0.5, legend=false, xlim=(0,1), ylim=(0,1), xlabel=L"x_{n}", ylabel=L"x_{n+1}")
end
mat_orb3
pog
savefig("lyapunov.pdf")
mixing = @animate for i in 2:n_orbits
    #scatter(1:size(mat_orb3, 1), mat_orb3[:, 1], label="Trajectory 1", color=:blue, markersize=1.5, markerstrokewidth =0.5, legend=false, xlim=(0,50), ylim=(0,1))
    scatter!(1:i,  mat_orb3[1:n, :], label="Trajectory $i", markersize=1.5, markerstrokewidth =0.5, legend=false, xlim=(0,50), ylim=(0,1))
end

mixing = @animate for n in 1:size(mat_orb3, 1)
    # Create a scatter plot for the current frame
    scatter(1:n, mat_orb3[1:n, :], 
            label=["Trajectory $i" for i in 1:n_cols],
            markersize=3, 
            xlabel="n", 
            ylabel="x_{n+1}", 
            title="Trajectories vs n", 
            legend=false, 
            xlim=(0, 50),
            ylim=(0,1),
            size=(600, 400))
end

gif(mixing, "mixing.gif", fps = 2)
