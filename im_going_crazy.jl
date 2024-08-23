include("methods/EI_estimate.jl")
include("methods/chaotic_system_methods.jl")


Random.seed!(1234)

n_orbits = 10
orbit_length = 7
initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)
window_sizes = collect(1:10)
a = 2
pertubation = 1/10^3
p0 = 0
p1 = 1/(sqrt(2))




testing = gumbel_window_min(mat_orb3, window_sizes)



gumby = observable_two(test, p0)
fit1 = gumbelfit(moving_minimum(gumby, 1))
fit2 = gumbelfit(moving_minimum(gumby, 2))
fit2 = gumbelfit(moving_minimum(gumby, 3))

p1 = location(fit1)
p2 = location(fit2) 