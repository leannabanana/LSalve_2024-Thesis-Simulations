include("chaotic_system_methods.jl")


import Plots as pl
using Plots, Extremes, Distributions, Gadfly, LaTeXStrings

set_default_plot_size(25cm, 20cm)
#CHAOS MAP
function chaotic_map(a, n_steps)
    y_values = Float64[]

    for i in 0:1/n_steps:1
        x = (a*i) % 1
        push!(y_values, x)

    end
    return y_values
end

function T_mod_map(start, n)
    y_values = Float64[]
    x = start

    for i in 0:n
        x = (3*x) % 1
        push!(y_values, x)

    end
    return y_values
end

function T_mod_map_noisey(start, n)
    y_values = Float64[]
    x = start
    normal_dist_t = Normal(0, 1/(10^6))

    for i in 0:n
        x = (3*x) % 1 + rand(normal_dist_t)
        push!(y_values, x)

    end
    return y_values
end

#Defining our observables
function observable_one(vector_values, x0, alpha)
    distances = abs.(vector_values .- x0).^(-alpha)
    return distances
end

function observable_two(vector_values, x0)
    distances = -log.(abs.(vector_values .- x0))
    return distances
end

function observable_three(vector_values, x0, a, alpha)
    distances = a.- (abs.(vector_values .- x0).^(-alpha))
    return distances
end


alpha = 1/6
x0 = 1/3
c = 3
start = 0.01
interations = 10^4

y = chaotic_map(3, interations)
x_values = range(0, 1, length=length(y))

y_vales_composed  = observable_one(T_mod_map_noisey(start, interations), x0 , alpha)
y_values_composed_2 = observable_two(T_mod_map_noisey(start, interations), x0 )
y_values_composed_3 = observable_three(T_mod_map_noisey(start, interations), x0 , c, alpha)

#### Generating Noise for our epic system
# Define the normal distribution with mean 0 and standard deviation 10^-6
pl.plot(x_values, y, xlimits=(0,1), legend=false)
p1 = scatter!(x_values, y_vales_composed, mc =:pink, ms=2, ma=0.5, title=L"\phi = d(x, x_0)^{-α}")
fittest=gevfit(y_vales_composed)
d111 = diagnosticplots(fittest)

pl.plot(x_values, y, xlimits=(0,1), legend=false)
p2 = scatter!(x_values, y_values_composed_2, mc =:pink, ms=2, ma=0.5, title=L"\phi = -log(d(x, x_0))")
fittest2 =gevfit(y_values_composed_2)
d222 = diagnosticplots(fittest2)

pl.plot(x_values, y, xlimits=(0,1), label = L"T(x) = 3x \, \, mod \, \, 1", mc=:purple, legend=false)
p3 = scatter!(x_values, y_values_composed_3, mc =:pink, ms=2, ma=0.5, title=L"\phi = c - d(x, x_0)^{-α}")
testtt = pl.plot(p1, p2, p3, layout=(1,3), size=(1100,400))
fittest3 =gevfit(y_values_composed_3)
d333 = diagnosticplots(fittest3)

iterate_y = y_vales_composed[1:end-1]
iterate_x =  y_vales_composed[2:end]
pl.plot(iterate_x, iterate_y)
scatter!(x_values, y_vales_composed, mc =:pink, ms=2, ma=0.5, title=L"\phi = d(x, x_0)^{-α}")

savefig(testtt,"newrealbigplot.png")
draw(SVG("obs1.svg", 25cm, 15cm), d111)

draw(SVG("obs2.svg",25cm, 15cm), d222)

draw(SVG("obs3.svg", 25cm, 15cm), d333)

function chaotic_map_test(x0, n_steps)
    x_values = Float64[]
    x = x0

    for i in 0:0.1:n_steps
        x = (3*x) % 1
        push!(x_values, x)
    end
    return x_values
end

y_vales_composed
x0 = 0.1
n_steps = 1

y_values1 = chaotic_map_test(x0, n_steps)
x_values1 = range(0, 1, length=length(x_values))

pl.plot(x_values1, y_values1)
scatter!(x_values, y_vales_composed, mc =:pink, ms=2, ma=0.5, title=L"\phi = d(x, x_0)^{-α}")

pl.plot(x_values, y, xlimits=(0,1), legend=false)
pl.plot()
pl.plot!(x_values, y_vales_composed, mc =:pink, ms=2, ma=0.5, title=L"\phi = d(x, x_0)^{-α}")


pl.plot(x_values, y, xlimits=(0,1), legend=false)
p2 = pl.plot!(x_values, y_values_composed_2, mc =:pink, ms=2, ma=0.5, title=L"\phi = -log(d(x, x_0))")



pl.scatter(x_values, y, xlimits=(0,1), label = L"T(x) = 3x \, \, mod \, \, 1", mc=:purple, legend=false)
p3 = pl.plot!(x_values, y_values_composed_3, mc =:pink, ms=2, ma=0.5, title=L"\phi = c - d(x, x_0)^{-α}")


using Plots

# Define the map function
T(x) = 3 * x % 1

# Number of iterations
N = 1000

# Range of initial conditions
x_range = collect(0:0.01:1)  # Adjust the step size for more or less granularity

# Arrays to store the data points
x_n = []
x_nplus1 = []

# Iterate the map for each initial condition
for x₀ in x_range
    x = x₀
    for i in 1:N
        x_next = T(x)
        push!(x_n, x)
        push!(x_nplus1, x_next)
        x = x_next
    end
end

# Plot x_{n+1} vs x_n
scatter(x_n, x_nplus1, markersize=2, xlabel="xₙ", ylabel="xₙ₊₁", legend=false)
pl.plot!(x_values, y)
title!("Graph of xₙ₊₁ vs xₙ for T(x) = 3x mod 1 map")
