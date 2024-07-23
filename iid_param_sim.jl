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
window_sizes = collect(1:13)
num_vectors = 10^3
vector_size = 10^4

distribution = Exponential(5)
X_n = generate_rv(num_vectors, vector_size, distribution)
X_n_mat = reduce(hcat, X_n)
X_n
### Define a function which gets me parameters for changing window sizes for the moving minimum functional 
function rv_minimum(random_variables, window_sizes)
    EI = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    Threads.@threads for windows in window_sizes
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

mus = scatter(window_sizes, iid_case[1], title=L"μ", legend=false, xlabel = "k", lc = "indianred2")
sigma = scatter(window_sizes, iid_case[2], title=L"σ", legend=false, xlabel = "k",  lc = "indianred2")
theta = scatter(window_sizes, iid_case[3], title=L"θ", legend=false, xlabel = "k",  lc = "indianred2")
min_params_1 = pl.plot(mus, sigma, theta, size=(900, 450), layout=(1,3), plot_title="Simluated RV iid rv vs window size moving minimum")

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
    return EI, location_params, scale_params
end

#More plots
iid_case_av = rv_average(X_n, window_sizes)

mus_av = pl.plot(window_sizes, iid_case_av[2],legend=false, ylabel=L"μ", xlabel=L"k")
sigma_av = pl.plot(window_sizes, iid_case_av[3], legend=false, ylabel=L"σ", xlabel=L"k")
theta_av = pl.plot(window_sizes, iid_case_av[1], legend=false, ylabel=L"σ", xlabel=L"k")
av_params = pl.plot(mus_av, sigma_av, theta_av, size=(900,500), layout=(1,3), plot_title="Simulated RV Moving Average Parameters vs Window Size")

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

    for i in 10:size(random_variables)[1]
        minimums = moving_minimum_matrix(random_variables[1:i, :], window_size)
        max_min = maximum(minimums, dims=1)[:]

        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
     return location_params, scale_params
end

what = test_3(X_n_mat, 10)
theme(:muted)


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

    Threads.@threads for i in 10:size(random_variables)[1]
        minimums = moving_average_matrix(random_variables[1:i, :], window_size)
        max_min = maximum(minimums, dims=1)[:]

        fit = gumbelfit(max_min)
        location_1 = location(fit)
        scale_1 = scale(fit)

        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
     return location_params, scale_params
end


test_4(X_n_mat, 10)


x_axis =  collect(10:size(X_n_mat)[1])
σ_k_av =scatter(x_axis, what_2[2], legend=false, xlabel = "Block Length", ylabel=L"\sigma^*", ms=1/2, ma =1/2, mc="indianred1", markerstrokecolor="indianred1", title ="GEV Maxima sampled from Exp(3)",  gridcolor=:gray19, gridalpha=1/2)
μ_k_av = scatter(x_axis, what_2[1], legend=false, xlabel = "Block Length", ylabel=L"\mu^*", ms=1/2, ma =1/2, mc="indianred1",  markerstrokecolor="indianred1",  title ="GEV Maxima sampled from Exp(3)", gridcolor=:gray19, gridalpha=1/2)

savefig(μ_k_av,"Output_Images/dependent_vs_independent_params/location_param_theorem_av.pdf")
savefig(σ_k_av,"Output_Images/dependent_vs_independent_params/scale_param_theorem_av.pdf")
