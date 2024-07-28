include("methods/chaotic_system_methods.jl")
include("methods/EI_estimate.jl")

x_values = collect(10:1:10^4)

av_iid = CSV.read("Data_csv/iid_av_blocklength.csv", DataFrame)
min_iid = CSV.read("Data_csv/iid_min_blocklength.csv", DataFrame)

av_reccurrent = CSV.read("Data_csv/av_length_vs_parameters_recurrent.csv", DataFrame)
av_nonrecurrent = CSV.read("Data_csv/av_length_vs_parameters_nonrecurrent.csv", DataFrame)

min_recurrent = CSV.read("Data_csv/min_length_vs_parameters_recurrent.csv", DataFrame)
min_nonrecurrent = CSV.read("Data_csv/min_length_vs_parameters_nonrecurrent.csv", DataFrame)
new_row = (location = 0, scale = 0, EI = 0)

# Insert the new row at the start (position 1) of the DataFrame
insert!(min_nonrecurrent, 1, new_row)


pl.plot(x_values, min_iid.location, label="Random Variable Case")
pl.plot!(x_values, min_nonrecurrent.location, label="nonrecurrent moving min")
pl.plot!(x_values, min_recurrent.location, label="recurrent moving min")

scatter(x_values, av_nonrecurrent.location,  markerstrokewidth = 1/9, ms=2/3)
pl.plot!(x_values, av_iid.location .+  (av_iid.scale ./ 2).*log.(av_nonrecurrent.EI ./ av_iid.EI))

pain =  ( exp.((av_nonrecurrent.location .- av_iid.location )./ av_nonrecurrent.scale) ).*av_nonrecurrent.EI

pl.plot(x_values, pain)
pl.plot!(x_values, av_iid.EI)

scatter(x_values, min_nonrecurrent.scale)




