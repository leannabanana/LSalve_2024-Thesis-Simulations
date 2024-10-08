include("chaotic_system_methods.jl")

### Estimate θ

function extremal_FerroSegers(Y, p)
    # This function computes the extremal index theta by using the
    # method proposed by Ferro-Segers (Ferro, C. A. T., and
    # J. Segers (2003), Inference for clusters of extremes,
    # J. R. Stat. Soc., Ser. B, 65, 545-556.).

    # INPUTS:
    # - Y: a vector containing a univariate time series
    # - p: a quantile value
    # OUTPUTS:
    # - theta: the estimate of the extremal index.
    
    # Extract the threshold u corresponding to the quantile p
    u = quantile(Y, p)
    
    # Compute the exceedances
    Si = findall(y -> y > u, Y)
    
    # Compute the cluster lengths
    Ti = diff(Si)
    
    # Compute the total number of clusters
    N = length(Ti)
    
    # Use the Ferro-Segers formula to extract theta
    theta = 2 * (sum(Ti .- 1))^2 / (N * sum((Ti .- 1) .* (Ti .- 2))) 
    
    return theta
end

# extremal_FerroSegers(gev_max, 0.95)

function EI_estimate_plot(orbits, window_size, x0)
    EI = Float64[]
    location = Float64[]
    scale = Float64[]

     for i in 2:length(orbits)
        observable = map(orbit -> observable_two(orbit, x0), orbits)
        min = moving_minimum.(observable, window_size)
        max = maximum.(min)

        EI_estimate = extremal_FerroSegers(max, 0.95)
        fit = gumbelfit(max)
        location_pm = location(fit)
        scale_pm = scale(fit)
    
        append!(location, location_pm)
        append!(scale, scale_pm)
        append!(EI, EI_estimate)
    end
    return EI, location, scale
end

# @btime EI_estimate_plot(orbits, 10, 0)
# println("x0 = 0 gets us 0.776479812607765 and x0 = 1/sqrt(2) = 1.0181")


function f_EI_window_min(orbits, window_sizes)
    location_params = Vector{Vector{Float64}}(undef, length(window_sizes))
    scale_params = Vector{Vector{Float64}}(undef, length(window_sizes))
    EI = Vector{Float64}(undef, length(window_sizes))
    shape_params = Vector{Vector{Float64}}(undef, length(window_sizes))
    
    Threads.@threads for i in 1:length(window_sizes)
        windows = window_sizes[i]
        minimums = moving_minimum_matrix(orbits, windows)
        max_min = maximum(minimums, dims=1)[:]

        fit = gumbelfit(max_min)

        EI_estimate = extremal_FerroSegers(max_min, 0.95)
        location_1 = location(fit)
        scale_1 = scale(fit)
        shapes = shape(fit)

        location_params[i] = location_1
        scale_params[i] = scale_1
        EI[i] = EI_estimate
        shape_params[i] = shapes
    end
    
    return location_params, scale_params, EI, shape_params
end

function EI_estimation_average(orbits, window_size)
    # location_params = Vector{Float64}(undef,19)
    # scale_params = Vector{Float64}(undef, 19)

    location_params = Float64[]
    scale_params = Float64[]
    EI = Float64[]

    for i in 10:size(orbits)[1]
        minimums = moving_average_matrix(orbits[1:i, :], window_size)
        max_min = maximum(minimums, dims=1)[:]

        EI_estimate = extremal_FerroSegers(max_min, 0.95)
        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)
        append!(EI, EI_estimate)
    end
     return location_params, scale_params, EI
end



function EI_estimation_min(orbits, window_size)
    # location_params = Vector{Float64}(undef,19)
    # scale_params = Vector{Float64}(undef, 19)

    location_params = Float64[]
    scale_params = Float64[]
    EI = Float64[]

 
    for i in 10:size(orbits)[1]
        minimums = moving_minimum_matrix(orbits[1:i, :], window_size)
        max_min = maximum(minimums, dims=1)[:]

        EI_estimate = extremal_FerroSegers(max_min, 0.95)
        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)
        append!(EI, EI_estimate)
    end
     return location_params, scale_params, EI
end

function EI_window_av(orbits, window_sizes)
    location_params = Float64[]
    scale_params = Float64[]
    EI = Float64[]
    shape_params = Float64[]

    
    for windows in window_sizes
        minimums = moving_average_matrix(orbits, windows)
        max_min = maximum(minimums, dims=1)[:]

        EI_estimate = extremal_FerroSegers(max_min, 0.95)
        location_1 = location(gevfit(max_min))
        scale_1 = scale(gevfit(max_min))
        shapes = shape(gevfit(max_min))

        append!(location_params, location_1)
        append!(scale_params, scale_1)
        append!(EI, EI_estimate)
        append!(shape_params, shapes)
    end
    return location_params, scale_params, EI, shape_params

end

function EI_window_min(orbits, window_sizes)
    location_params = Float64[]
    scale_params = Float64[]
    EI = Float64[]
    shape_params = Float64[]
    
    for windows in window_sizes
        minimums = moving_minimum_matrix(orbits, windows)
        max_min = maximum(minimums, dims=1)[:]

        fit = gevfit(max_min)

        EI_estimate = extremal_FerroSegers(max_min, 0.95)
        location_1 = location(fit)
        scale_1 = scale(fit)
        shapes = shape(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)
        append!(EI, EI_estimate)
        append!(shape_params, shapes)
    end
    return location_params, scale_params, EI, shape_params

end


function frequency_plots_min(data, window_sizes)
    for windows in window_sizes

        minimums = moving_minimum_matrix(data, windows)
        max_min = maximum(minimums, dims=1)[:]

        fit = gevfit(max_min)
        freq = histplot(fit)
        draw(PDF("Output_Images/frequency_plots/frequency_window"*string(windows)*".pdf",25cm, 15cm), freq)

    end
end



function gumbel_window_min(orbits, window_sizes)
    location_params = Float64[]
    scale_params = Float64[]
    EI = Float64[]

    for windows in window_sizes
        minimums = moving_minimum_matrix(orbits, windows)
        max_min = maximum(minimums, dims=1)[:]
    
        EI_estimate = extremal_FerroSegers(max_min, 0.95)
        location_1 = location(gumbelfit(max_min))
        scale_1 = scale(gumbelfit(max_min))
  
    
        append!(location_params, location_1)
        append!(scale_params, scale_1)
        append!(EI, EI_estimate)
    
        end
        return location_params, scale_params, EI
end


function gumbel_window_av(orbits, window_sizes)
    location_params = Float64[]
    scale_params = Float64[]
    EI = Float64[]

    for windows in window_sizes
        minimums = moving_average_matrix(orbits, windows)
        max_min = maximum(minimums, dims=1)[:]
    
        EI_estimate = extremal_FerroSegers(max_min, 0.95)
        location_1 = location(gumbelfit(max_min))
        scale_1 = scale(gumbelfit(max_min))
    
        append!(location_params, location_1)
        append!(scale_params, scale_1)
        append!(EI, EI_estimate)
    
        end
        return location_params, scale_params, EI
end


function leanna_mu( λ, k, μ_1, σ_1, θ_2, θ_1) ### returns μ_2

    μ_2 = λ^(-(k-1.0))*μ_1 + λ^(-(k-1.0))*σ_1*log(θ_2/ θ_1)

    return μ_2
end

function leanna_mu_2( k, μ_1, σ_1, θ_2, θ_1 )
    μ_2 = μ_1/k + ((σ_1 /k) *log(θ_2/ θ_1)) 

    return μ_2
end
