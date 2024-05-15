"""
This file verifies our number of blocks
"""

include("GEV_fitting.jl")
include("chaotic_system_methods.jl")

k = collect(9:150) 

shapes_1 = verify_blocks(y_obs_1)
xi_1 = pl.plot(k, shapes_1, legend=false, title=L"\phi = d(x, x_0)^{-α}")

shapes_2 = verify_blocks(y_obs_2)
xi_2 = pl.plot(k, shapes_2, legend=false, title=L"\phi = -log(d(x, x_0))")

shapes_3 = verify_blocks(y_obs_3)
xi_3 = pl.plot(k, shapes_3, legend=false, title=L"\phi = c - d(x, x_0)^{-α}")

k_plots = pl.plot(xi_1, xi_2, xi_3, layout=(1,3), ylabel=L"ξ",  xlabel="Number of k blocks", size=(1300,500))
savefig(k_plots,"Output_Images/gev_diagnostic_tests/blocks_vs_xi.png")
