"""
This file puts all the outputs into a df and exports it as an excel document
"""

include("GEV_fitting.jl")
include("methods/chaotic_system_methods.jl")

# Create DataFrame
observables = DataFrame(Perturbed_Map = y_n_noise, Observable_1_Frechet = y_obs_1, Observable_2_Gumbel = y_obs_2, Observable_3_Weibull= y_obs_3)

maximums = DataFrame(Max_Frechet = maximums_1, Max_Gumbel = maximums_2, Max_Weibul = maximums_3)

CSV.write("Data_csv/observables_data.csv", observables, delim=',', header=true)
CSV.write("Data_csv/maximum_observables_data.csv", maximums, delim=',', header=true)