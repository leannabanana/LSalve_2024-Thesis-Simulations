"""
This file outputs simulations + fits GEV according to methods in "chaotic_system_methods.jl"
"""

include("chaotic_system_methods.jl")
Date = "_07-05-24_10^-6"
Random.seed!(1234)
#Define Constants
a = 3
alpha = 1/6
x0 = 1/3
c = 1
interations = 10^4
pertubation = 1/(10^6)

#Create T(x) = 3xmod1
y_map = chaotic_map(a, interations)
x_map = range(0, 1, length=length(y_map))
pl.plot(x_map, y_map)

### T^n(x)
y_n = T_mod_map_noisey(a, interations, pertubation)

### Perturbed Chaotic Map 3x mod 1 with our 3 observables
y_obs_1 = observable_one(T_mod_map_noisey(a, interations, pertubation), x0 , alpha)
y_obs_2 = observable_two(T_mod_map_noisey(a, interations, pertubation), x0 )
y_obs_3 = observable_three(T_mod_map_noisey(a, interations, pertubation), x0 , c, alpha)

### We plot the scatter plots onto T(x)
#Observable 1
 
scatter(y_n[1:end-1],y_n[2:end], ms=2, ma=0.5, mc=:lightblue,  xlims=(0,1), legend=false)
p1 = scatter!(x_map, y_obs_1, mc =:pink, ms=2, ma=0.5, title=L"\phi = d(x, x_0)^{-α}")

#Observable 2
scatter(y_n[1:end-1],y_n[2:end], ms=2, ma=0.5, mc=:lightblue, xlims=(0,1), legend=false)
p2 = scatter!(x_map, y_obs_2, mc =:pink, ms=2, ma=0.5, title=L"\phi = -log(d(x, x_0))")

#Observable 3
scatter(y_n[1:end-1],y_n[2:end], ms=2, ma=0.5, mc=:lightblue, xlims=(0,1), legend=false)
p3 = scatter!(x_map, y_obs_3, mc =:pink, ms=2, ma=0.5, title=L"\phi = c - d(x, x_0)^{-α}")

### Put 3 plots onto one image 
big_plot = pl.plot(p1, p2, p3, layout=(1,3), size=(1100,400))
savefig(big_plot,"Output_Images/observable_scatter_plots/3_observables_10^-6_perturbation.png")

### Fitting a GEV
#Observable 1
fit_obs1 = gevfit(y_obs_1)
d1 = diagnosticplots(fit_obs1)

#Observable 2
fit_obs2 =gevfit(y_obs_2)
d2 = diagnosticplots(fit_obs2)

#Observable 3
fit_obs3 =gevfit(y_obs_3)
d3 = diagnosticplots(fit_obs3)

#Save Diagnostic Tests
draw(PDF("Output_Images/gev_diagnostic_tests/obs1"*Date*".pdf", 25cm, 15cm), d1)
draw(PDF("Output_Images/gev_diagnostic_tests/obs2"*Date*".pdf",25cm, 15cm), d2)
draw(PDF("Output_Images/gev_diagnostic_tests/obs3"*Date*".pdf", 25cm, 15cm), d3)