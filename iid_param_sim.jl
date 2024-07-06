include("chaotic_system_methods.jl")

Random.seed!(1234)

### Define window sizes
window_sizes = collect(1:50)

### Simulate a bunch of RV's 
function exponential_distributions(n_vectors, size, distribution)
    random_vectors = [rand(distribution, vector_size) for _ in 1:num_vectors]
    return random_vectors
end

num_vectors = 10^3
vector_size = 10^3
λ = 12 #just because i can
distribution = Exponential(λ)

X_n = exponential_distributions(num_vectors, vector_size, distribution)

### Define a function which gets me parameters for changing window sizes for the moving minimum functional 
function rv_minimum(random_variables, window_sizes)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    for windows in window_sizes
        minimums = moving_minimum.(random_variables, windows)
        max_min = maximum.(minimums)
        
        shapes = 0.0
        location_1 = 0.0
        scale_1 = 0.0
        
        try
            fit = gumbelfit(max_min)
            shapes = shape(fit)
            location_1 = location(fit)
            scale_1 = scale(fit)
        catch e
            if isa(e, DomainError)
                shapes = 0.0
                location_1 = 0.0
                scale_1 = 0.0
            else
                rethrow(e)
            end
        end

        append!(shape_params, shapes)
        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
    return shape_params, location_params, scale_params
end

#plots
iid_case = rv_minimum(X_n, window_sizes)
mus = pl.plot(window_sizes, iid_case[2], title=L"μ", legend=false)
theta = pl.plot(window_sizes, iid_case[3], title=L"θ", legend=false)
min_params_1 = pl.plot(mus, theta, size=(800, 600), layout=(1,2), plot_title="Simluated RV iid rv vs window size moving minimum")


### Define a function which gets me parameters for changing window sizes for the moving average functional
function rv_average(random_variables, window_sizes)
    shape_params = Float64[]
    location_params = Float64[]
    scale_params = Float64[]
    
    for windows in window_sizes
        minimums = moving_average.(random_variables, windows)
        max_min = maximum.(minimums)
        
        shapes = 0.0
        location_1 = 0.0
        scale_1 = 0.0
        
        try
            fit = gumbelfit(max_min)
            shapes = shape(fit)
            location_1 = location(fit)
            scale_1 = scale(fit)
        catch e
            if isa(e, DomainError)
                shapes = 0.0
                location_1 = 0.0
                scale_1 = 0.0
            else
                rethrow(e)
            end
        end

        append!(shape_params, shapes)
        append!(location_params, location_1)
        append!(scale_params, scale_1)
    end
    return shape_params, location_params, scale_params
end

#More plots
iid_case_av = rv_average(X_n, window_sizes)
mus_av = pl.plot(window_sizes, iid_case_av[2],legend=false, ylabel=L"μ", xlabel=L"k")
theta_av = pl.plot(window_sizes, iid_case_av[3], legend=false, ylabel=L"σ", xlabel=L"k")


av_params = pl.plot(mus_av, theta_av, size=(800,600), layout=(1,2), plot_title="Simulated RV Moving Average Parameters vs Window Size")
av_params_min = pl.plot(iid_case[2] .- iid_case[3].*log.(window_sizes), iid_case[2], legend=false, ylabel=L"\mu_k", xlabel=L"\mu = \mu_k - σ\log(k)", title=L"$\mu_k$ vs $\mu$ (moving minimum)")
av_params_av = pl.plot(iid_case_av[2] .- iid_case_av[3].*log.(window_sizes), iid_case_av[2], legend=false, ylabel=L"\mu_k", xlabel=L"\mu = \mu_k - σ\log(k)", title=L"$\mu_k$ vs $\mu$  (moving average)")

#Saving the plots
savefig(av_params,"Output_Images/verifying_dependent_iid_parameters/Moving_average_params_vs_windowsize.pdf")
savefig(min_params_1,"Output_Images/verifying_dependent_iid_parameters/Moving_min_params_vs_windowsize.pdf")
savefig(av_params_min,"Output_Images/verifying_dependent_iid_parameters/muk_vs_mu_movingmin.pdf")
savefig(av_params_av,"Output_Images/verifying_dependent_iid_parameters/muk_vs_mu_movingav.pdf")