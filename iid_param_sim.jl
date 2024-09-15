"""
This file simulates iid random variables to check how Gumbel parameters change
"""

include("methods/chaotic_system_methods.jl")
include("methods/EI_estimate.jl")

Random.seed!(1234)

### Simulate a bunch of RV's 
function generate_rv(n_vectors, vector_size, distribution)
    random_vectors = [rand(distribution, vector_size) for _ in 1:n_vectors]
    return random_vectors
end


### Define window sizes
window_sizes = collect(1:10)
num_vectors = 10^3
vector_size = 10^3

distribution = Exponential(5)

X_n = generate_rv(num_vectors, vector_size, distribution)
X_n_mat = reduce(hcat, X_n)

### Define a function which gets me parameters for changing window sizes for the moving minimum functional 
function rv_minimum(random_variables, window_sizes)
    EI = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    for windows in window_sizes
        minimums = moving_minimum.(random_variables, windows)
        max_min = maximum.(minimums)
        EI_estimate = extremal_FerroSegers(max_min, 0.95)
        fit = gumbelfit(max_min)
        #shapes = shape(fit)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(EI, EI_estimate)
        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
    return  location_params, scale_params, EI
end

#plots
iid_case = rv_minimum(X_n, window_sizes)

mus = scatter(window_sizes, iid_case[1], title="Gumbel Distribution Sampled from Exp(5)", legend=false, ylabel =L"\mu", xlabel = "k", mc = "indianred2", ms = 2.5)
pl.plot!(window_sizes, iid_case[1][1] ./ window_sizes .+ iid_case[2][1].*log.(iid_case[3][1] ./ iid_case[3]) )

sigma = scatter(window_sizes, iid_case[2], title="Gumbel Distribution Sampled from Exp(5)", legend=false, xlabel = "k", ylabel=L"\sigma",  mc = "indianred2", ms = 2.5)
pl.plot!(window_sizes, iid_case[2][1]./window_sizes)


savefig(mus, "Output_Images/finalised_iid_plots/minimum_mu.pdf")
savefig(sigma, "Output_Images/finalised_iid_plots/minimum_sigma.pdf")

theta = scatter(window_sizes, iid_case[3], title=L"θ", legend=false, xlabel = "k",  mc = "indianred2", ms = 2.5)
min_params_1 = pl.plot(mus, sigma, size=(650, 400), layout=(1,2), plot_title="Simluated Exp(3) vs window size moving minimum")
savefig(min_params_1,"Output_Images/verifying_dependent_iid_parameters/aaaaa.pdf")

### Define a function which gets me parameters for changing window sizes for the moving average functional
function rv_average(random_variables, window_sizes)
    EI = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    for windows in window_sizes
        minimums = moving_average.(random_variables, windows)
        max_min = maximum.(minimums)
        EI_estimate_plot = extremal_FerroSegers(max_min, 0.95)

        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)
        append!(EI, EI_estimate_plot)
        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
    return location_params, scale_params, EI
end

#More plots
iid_case_av = rv_average(X_n, window_sizes)
scatter(window_sizes, iid_case_av[1])
pl.plot!(window_sizes, 1.5 .* iid_case_av[1][1] ./( window_sizes .+ 1 ))


mu_av = scatter(window_sizes, iid_case_av[1], title="Gumbel GEV Sampled from Exp(5) - Moving Average", legend=false, ylabel =L"\mu", xlabel = "k", mc = "indianred2", ms = 2.5)
pl.plot!(window_sizes, iid_case_av[1][1]./(window_sizes) )


sigma_av = scatter(window_sizes, iid_case_av[2],  title="Gumbel GEV Sampled from Exp(5) - Moving Average", legend=false, ylabel =L"\sigma", xlabel = "k", mc = "indianred2", ms = 2.5)
pl.plot!(window_sizes, iid_case_av[2][1] ./ window_sizes )

savefig(mu_av, "Output_Images/finalised_iid_plots/av_mu.pdf")
savefig(sigma_av, "Output_Images/finalised_iid_plots/minimum_sigma.pdf")



theta_av = pl.plot(window_sizes, iid_case_av[1], legend=false, ylabel=L"σ", xlabel=L"k")
av_params = pl.plot(mus_av, sigma_av, size=(650, 400), layout=(1,2), plot_title="Simluated Exp(3) vs window size moving average")

av_params_min = pl.plot(iid_case[2] .- iid_case[3][1].*log.(window_sizes), iid_case[2], legend=false, ylabel=L"\mu_k", xlabel=L"\mu = \mu_k - σ\log(k)", title=L"$\mu_k$ vs $\mu$ (moving minimum)")
av_params_av = pl.plot(iid_case_av[2] .- iid_case_av[3][1].*log.(window_sizes), iid_case_av[2], legend=false, ylabel=L"\mu_k", xlabel=L"\mu = \mu_k - σ\log(k)", title=L"$\mu_k$ vs $\mu$  (moving average)")
av_params_min_k = pl.plot(iid_case[2] .- iid_case[3].*log.(window_sizes), iid_case[2], legend=false, ylabel=L"\mu_k", xlabel=L"\mu = \mu_k - σ_k\log(k)", title=L"$\mu_k$ vs $\mu$ (moving minimum)")
av_params_av_k  = pl.plot(iid_case_av[2] .- iid_case_av[3].*log.(window_sizes), iid_case_av[2], legend=false, ylabel=L"\mu_k", xlabel=L"\mu = \mu_k - σ_k\log(k)", title=L"$\mu_k$ vs $\mu$  (moving average)")
combined = pl.plot(av_params_min, av_params_av, av_params_min_k, av_params_av_k, layout = (2,2), size=(700, 550))

#Saving the plots
savefig(av_params,"Output_Images/verifying_dependent_iid_parameters/Moving_average_params_vs_windowsize.pdf")
savefig(min_params_1,"Output_Images/verifying_dependent_iid_parameters/Moving_min_params_vs_windowsize.pdf")
savefig(av_params_min,"Output_Images/verifying_dependent_iid_parameters/muk_vs_mu_movingmin.pdf")
savefig(av_params_av,"Output_Images/verifying_dependent_iid_parameters/muk_vs_mu_movingav.pdf")
savefig(combined,"Output_Images/verifying_dependent_iid_parameters/combined_plots.pdf")
3*10^4-1
X_n[1]

function test(random_variables, window_size)
    # location_params = Vector{Float64}(undef,19)
    # scale_params = Vector{Float64}(undef, 19)

    location_params = Float64[]
    scale_params = Float64[]

    Threads.@threads for i in 10:length(random_variables[1])
        data = [variable[1:i] for variable in random_variables]
        minimums = @. moving_minimum(data, window_size)
        max_min = @. maximum(minimums)

        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)

    end
     return location_params, scale_params
end

@btime this_sucks = test(X_n, 5)

pl.plot(aaa, testing[1])
pl.plot(aaa, testing[2])


max_values = maximum(min, dims=1)[:]

function test_2(random_variables, window_size)
    # location_params = Vector{Float64}(undef, 100)
    # scale_params = Vector{Float64}(undef, 100)

    location_params = Float64[]
    scale_params = Float64[]

    Threads.@threads for i in 2:length(random_variables)
        minimums = moving_minimum.(random_variables[1:i], window_size)
        max_min = maximum.(minimums)

        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)
    
    end
    return location_params, scale_params
end


function test_3(random_variables, window_size)
    # location_params = Vector{Float64}(undef,19)
    # scale_params = Vector{Float64}(undef, 19)

    location_params = Float64[]
    scale_params = Float64[]
    EI = Float64[]
    for i in 10:size(random_variables)[1]
        minimums = moving_minimum_matrix(random_variables[1:i, :], window_size)
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

x_axis =  collect(10:size(X_n_mat)[1])
σ_k =scatter(x_axis, what[2], legend=false, xlabel = "Block Length", ylabel=L"\sigma^*", ms=1/2, ma =1/2, mc="indianred1", markerstrokecolor="indianred1", title ="GEV Maxima sampled from Exp(3)",  gridcolor=:gray19, gridalpha=1/2)
μ_k = scatter(x_axis, what[1], legend=false, xlabel = "Block Length", ylabel=L"\mu^*", ms=1/2, ma =1/2, mc="indianred1",  markerstrokecolor="indianred1",  title ="GEV Maxima sampled from Exp(3)", gridcolor=:gray19, gridalpha=1/2)

savefig(μ_k,"Output_Images/dependent_vs_independent_params/location_param_theorem.pdf")
savefig(σ_k,"Output_Images/dependent_vs_independent_params/scale_param_theorem.pdf")




function test_4(random_variables, window_size)
    # location_params = Vector{Float64}(undef,19)
    # scale_params = Vector{Float64}(undef, 19)

    location_params = Float64[]
    scale_params = Float64[]
    EI = Float64[]

    for i in 10:size(random_variables)[1]
        minimums = moving_average_matrix(random_variables[1:i, :], window_size)
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


what = test_3(X_n_mat, 10)

what1 = DataFrame(location = what[1], scale =  what[2], EI =  what[3])
CSV.write("Data_csv/iid_min_blocklength.csv", what1, delim=',', header=true)


what2 = test_4(X_n_mat, 2)
what12 = DataFrame(location = what2[1], scale =  what2[2], EI =  what2[3])
CSV.write("Data_csv/iid_av_blocklength.csv", what12, delim=',', header=true)


x_axis =  collect(10:size(X_n_mat)[1])
σ_k_av =scatter(x_axis, what_2[2], legend=false, xlabel = "Block Length", ylabel=L"\sigma^*", ms=1/2, ma =1/2, mc="indianred1", markerstrokecolor="indianred1", title ="GEV Maxima sampled from Exp(3)",  gridcolor=:gray19, gridalpha=1/2)
μ_k_av = scatter(x_axis, what_2[1], legend=false, xlabel = "Block Length", ylabel=L"\mu^*", ms=1/2, ma =1/2, mc="indianred1",  markerstrokecolor="indianred1",  title ="GEV Maxima sampled from Exp(3)", gridcolor=:gray19, gridalpha=1/2)

savefig(μ_k_av,"Output_Images/dependent_vs_independent_params/location_param_theorem_av.pdf")
savefig(σ_k_av,"Output_Images/dependent_vs_independent_params/scale_param_theorem_av.pdf")

gumbel_dist = Gumbel(0, 1)

# Dimensions of the matrix
rows = 10^3
cols = 10^4

# Simulate the matrix
matrix_gumbel = rand(gumbel_dist, rows, cols)


iid = EI_window_av(matrix_gumbel, window_sizes)
scatter(window_sizes, iid[1],  xticks=1:1:13, xlabel = " k ", ylabel = L"μ", title=L"Moving minimum $x_0 = 0$", mc="tomato2", ms=3, ma=1)
pl.plot!(window_sizes, 2 .^log.(iid[1][1]./ window_sizes) .+ 2)


iid2 = EI_window_min(matrix_gumbel, window_sizes)
scatter(window_sizes, iid2[1],  xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title="Moving Minimum from an iid Gumbel distribution", mc="tomato2", ms=3, ma=1, legend=false)
pl.plot!(window_sizes, iid2[1][1] ./ window_sizes .+ iid2[2][1].*log.(iid2[3] ./ iid2[3][1]))


scatter(window_sizes, iid2[2],  xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title="Moving Minimum from an iid Gumbel distribution", mc="tomato2", ms=3, ma=1, legend=false)
pl.plot!(window_sizes, iid2[2][1] ./ window_sizes )

