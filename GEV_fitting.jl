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
initial_value = 0.01

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

### Plit into k_blocks and find its maximum_values
maximums_1 = maximum_values(y_obs_1, 50)
maximums_2 = maximum_values(y_obs_2, 50)
maximums_3 = maximum_values(y_obs_3, 50)

### We plot the scatter plots onto T(x)
#Observable 1
orbit_1 = scatter(y_n_noise[1:end-1],y_n_noise[2:end], ms=2, ma=0.5, mc=:lightblue,  xlims=(0,1), legend=false)
p1 = scatter!(x_map, y_obs_1, mc =:pink, ms=2, ma=0.5, title=L"\phi = d(x, x_0)^{-α}")

#Observable 2
orbit_2 = scatter(y_n_noise[1:end-1],y_n_noise[2:end], ms=2, ma=0.5, mc=:lightblue, xlims=(0,1), legend=false)
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



using Distributions
using Plots
using Late
# Define the parameters for each distribution
μ_gumbel, β_gumbel = 0.0, 1.0  # Gumbel distribution parameters
λ_weibull, k_weibull = 1.0, 1.5  # Weibull distribution parameters
α_frechet, s_frechet = 2.0, 1.0  # Fréchet distribution parameters

# Create the distributions
gumbel_dist = Gumbel(μ_gumbel, β_gumbel)
weibull_dist = Weibull(k_weibull, λ_weibull)
frechet_dist = Frechet(α_frechet, s_frechet)
GeneralizedExtremeValue(0, 1, 1)  



x = -5:0.01:5  # Ensure x > 0 for Weibull and Fréchet distributions
pdf_frechet = pdf.(GeneralizedExtremeValue(0, 1, 1/2)  , x)
pdf_gumbel = pdf.(GeneralizedExtremeValue(0, 1, 0)  , x)
pdf_weibull = pdf.(GeneralizedExtremeValue(0, 1, -1/2)  , x)

# Plot the PDFs
hi = plot(x, pdf_frechet, label="Gumbel", linewidth=2, lc="gold1", legend=:topright, xlabel=L"x", ylabel=L"f(x)")
plot!(x, pdf_gumbel, label="Weibull", linewidth=2, lc="lightslateblue")
plot!(x, pdf_weibull, label="Fréchet", linewidth=2, lc="palevioletred1")

title!("PDFs of Gumbel, Weibull, and Fréchet Distributions")
savefig("gumbel_weibull_frechet_pdfs.pdf")