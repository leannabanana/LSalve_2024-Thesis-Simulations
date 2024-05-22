"""
This file verifies our number of blocks
"""

include("GEV_fitting.jl")
include("chaotic_system_methods.jl")
include("GEV_fitting_functionals.jl")

k_block_number = collect(9:length(maximums_1)) 

### Verify Normal GEV fits
shapes_1 = xi_params(maximums_1, k_block_number)
xi_1 = pl.plot(k_block_number, shapes_1, legend=false, title=L"\phi = d(x, x_0)^{-α}")

shapes_2 = xi_params(maximums_2, k_block_number)
xi_2 = pl.plot(k_block_number, shapes_2, legend=false, title=L"\phi = -log(d(x, x_0))")

shapes_3 = xi_params(maximums_3, k_block_number)
xi_3 = pl.plot(k_block_number, shapes_3, legend=false, title=L"\phi = c - d(x, x_0)^{-α}")

k_plots = pl.plot(xi_1, xi_2, xi_3, layout=(1,3), ylabel=L"ξ",  xlabel="Number of k blocks", size=(1300,500))
savefig(k_plots,"Output_Images/gev_diagnostic_tests/blocks_vs_xi.png")


### Verify Functionals
fun_block_number = collect(4: length(av_max_1))
fun_1 = verify_blocks(av_max_1, fun_block_number)
fun_xi = pl.plot(fun_block_number, fun_1, legend=false, title=L"\phi = d(x, x_0)^{-α}")

fun_2 = verify_blocks(av_max_2, fun_block_number)
fun_xi_2 = pl.plot(fun_block_number, fun_2, legend=false, title=L"\phi = d(x, x_0)^{-α}")

fun_3 = verify_blocks(av_max_3, fun_block_number)
fun_xi_3 = pl.plot(fun_block_number, fun_3, legend=false, title=L"\phi = d(x, x_0)^{-α}")
