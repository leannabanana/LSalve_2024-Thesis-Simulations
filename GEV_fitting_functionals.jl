""""
This file Implements the moving average and moving minimum to the data from GEV_fitting.jl
"""

include("chaotic_system_methods.jl")
block_size = 300

## Moving Average
moving_av_1 = moving_average(y_obs_1, block_size)
moving_av_2 = moving_average(y_obs_2, block_size)
moving_av_3 = moving_average(y_obs_3, block_size)

av_max_1 = maximum_values(moving_av_1, 15)
av_max_2 = maximum_values(moving_av_2, 15)
av_max_3 = maximum_values(moving_av_3, 15)

gevfit(av_max_1)
gevfit(av_max_2)
gevfit(av_max_3)

### Moving Minimum
moving_min_1 = moving_minimum(y_obs_1, block_size)
moving_min_2 = moving_minimum(y_obs_2, block_size)
moving_min_3 = moving_minimum(y_obs_3, block_size)

min_max_1 = maximum_values(moving_av_1, 15)
min_max_2 = maximum_values(moving_av_2, 15)
min_max_3 = maximum_values(moving_av_3, 15)

gevfit(min_max_1)
gevfit(min_max_2)
gevfit(min_max_3)
