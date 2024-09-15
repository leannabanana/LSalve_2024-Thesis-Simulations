include("methods/EI_estimate.jl")
include("methods/chaotic_system_methods.jl")

Random.seed!(1234)

n_orbits = 10^3
orbit_length =  10^5
initial_conditions = collect(1/n_orbits : 1/n_orbits : 1)
window_sizes = collect(1:10)
a = 2
pertubation = 1/10^3
p0 = 0
p1 = 1/(sqrt(2))

orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
average_orbits = observable_two.(orbits, p1)
mat_orb3 = reduce(hcat, average_orbits)
testing = EI_window_av(mat_orb3, window_sizes)


df10 = DataFrame(location = testing[1], scale = testing[2], EI = testing[3])




simulate_orbits(0.001, 2, 10^6, pertubation, 1)

df10 = CSV.read("Data_csv/updated_data/final_gevfit_min_0.csv", DataFrame)

scatter(window_sizes, df10.location)
pl.plot!(window_sizes, log.(df10.location[1] .+ log.(2 .^(   (window_sizes .-1) )) .* ((df10.EI ./ df10.EI[1]).^(df10.scale[1]))   )   )


scatter(window_sizes, testing[2])
pl.plot!(window_sizes, testing[2][1] ./ window_sizes)

testing[2]


scatter(window_sizes, est[2])
pl.plot!(window_sizes, est[2][1]  ./ window_sizes)


df4 = DataFrame(location = av_est1[1], scale = av_est1[2], EI = av_est1[3], shape = av_est[4])
df4 = CSV.read("Data_csv/updated_data/av_0.csv", DataFrame)
scatter(window_sizes, df4.scale)
pl.plot!(window_sizes, df4.scale[1] ./ window_sizes)

scatter(window_sizes, log.(df4.location))
plot(window_sizes, df4.location[1])

df4.location[1]

estimate = EI_window_min(mat_orb, window_sizes)

average_orbits_1= observable_two.(orbits, 0)
mat_orb_1 = reduce(hcat, average_orbits_0)
estimate_0 = EI_window_min(mat_orb_1, window_sizes)


dfs = CSV.read("Data_csv/updated_data/final_gevfit_min_0.csv", DataFrame)
df1 = CSV.read("Data_csv/updated_data/final_gevfit_min_1sqrt2_data.csv", DataFrame)

new_mu = leanna_mu_2.(window_sizes, estimate[1][1], estimate[2][1], estimate[4], estimate[4][1])

dfmatrix = DataFrame(mat_orb, :auto)
CSV.write("Data_csv/updated_data/trajectories_106_1sqrt2.csv", estimate)

df = DataFrame(location = estimate[1], scale = estimate[2], EI = estimate[3], shape = estimate[4])
CSV.write("Data_csv/updated_data/final_gevfit_min_1sqrt2_data.csv", df, delim=',', header=true)

df = DataFrame(location = estimate_0[1], scale = estimate_0[2], EI = estimate_0[3], shape = estimate_0[4])
CSV.write("Data_csv/updated_data/final_gevfit_min_0.csv", df, delim=',', header=true)

mus = 2 .^( .- (window_sizes .-1.0)).*dfs.location[1] + dfs.scale[1]* 2 .^( .- (window_sizes .-1.0)).*log.(dfs.EI ./ dfs.EI[1])


g1 = scatter(window_sizes, dfs.location, xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title=L"Moving minimum $x_0 = 0$", mc="tomato2",  ms=3, ma=1)
pl.plot!(window_sizes, dfs.location[1] ./ 2 .^(window_sizes .- 1.0) + dfs.scale[1] .* log.(dfs.EI ./ dfs.EI[1]) ./ 2 .^(window_sizes .- 1.0))


g2 = scatter(window_sizes, dfs.scale, xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title=L"Moving minimum $x_0 = 0$", mc="tomato2", ms=3, ma=1)
pl.plot!(window_sizes, dfs.scale[1] ./ 2 .^(  (window_sizes .-1)  ))


g3 = scatter(window_sizes, dfs.EI, xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title=L"Moving minimum $x_0 = 0$", mc="tomato2", legend=false, ms=3, ma=1)
g4 = scatter(window_sizes, dfs.shape, xticks=1:1:13, xlabel = " k ", ylabel = L"ξ", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)

min = scatter(window_sizes, df1.location, xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title="Moving minimum non-recurrent", mc="tomato2",  ms=3, ma=1)
pl.plot!(window_sizes, df1.location[1] ./ 3 .^(window_sizes .- 1.0) + dfs.scale[1] .* log.(dfs.EI ./ dfs.EI[1]) ./ 3 .^(window_sizes .- 1.0))


min2 = scatter(window_sizes, df1.scale, xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title="Moving minimum non-recurrent", mc="tomato2", ms=3, ma=1)
pl.plot!(window_sizes, dfs.scale[1] ./ exp(2) .^(window_sizes .-1 ))


min3 = scatter(window_sizes, dfs.shape, xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title="Moving minimum non-recurrent", mc="tomato2", legend=false, ms=3, ma=1)
min4 = scatter(window_sizes, estimate[4], xticks=1:1:13, xlabel = " k ", ylabel = L"ξ", title="Moving minimum non-recurrent", mc="tomato2", legend=false, ms=3, ma=1)



orbits = simulate_orbits(initial_conditions, a, orbit_length, pertubation, n_orbits)
average_orbits_0 = observable_two.(orbits, 0)


mat_orb = reduce(hcat, average_orbits)
gumbelestimate = gumbel_window_min(mat_orb, window_sizes)
e1 = scatter(window_sizes, gumbelestimate[1], xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
pl.plot!(window_sizes, estimate[1][1].*(1 ./ 2 .^(window_sizes.-1.0)))

e2 = scatter(window_sizes, gumbelestimate[2], xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
e3 = scatter(window_sizes, gumbelestimate[3], xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title=L"Moving minimum $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)

gumbelmin = pl.plot(e1, e2, e3, layout = @layout([A B ; C D]), size=(700,550), plot_title="Gumbel Fit", plot_titlevspan=0.04)
savefig(gumbelmin,"Output_Images/18-08-2024/gumbel_min_1sqrt2.pdf")

m1 = scatter(window_sizes, estimate_av[1], xticks=1:1:13, xlabel = " k ", ylabel = L"\mu", title=L"Moving av $x_0 = 0$", mc="tomato2", legend=false, ms=3, ma=1)
pl.plot!(window_sizes, estimate_av[1][1] ./ window_sizes .+ (38.2271057231807 ./window_sizes).*log.(estimate_av[3] ./ estimate_av[3][1]))




m2 = scatter(window_sizes, estimate_av[2], xticks=1:1:13, xlabel = " k ", ylabel = L"\sigma", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
pl.plot!(window_sizes, estimate_av[2][1] ./ window_sizes)

m3 = scatter(window_sizes, estimate_av[3], xticks=1:1:13, xlabel = " k ", ylabel = L"\theta", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)
m4 = scatter(window_sizes, estimate_av[4], xticks=1:1:13, xlabel = " k ", ylabel = L"ξ", title=L"Moving av $x_0 = \dfrac{1}{\sqrt{2}}$", mc="tomato2", legend=false, ms=3, ma=1)

avfit = pl.plot(m1, m2, m3, m4, layout = @layout([A B ; C D]), size=(700,550), plot_title="GEV Fit", plot_titlevspan=0.04)
savefig(avfit,"Output_Images/updated_plots/gevfit_av__1sqrt2.pdf")

################ USING ONLY CSV 

m_0 = CSV.read("Data_csv/updated_data/min_0_fix.csv", DataFrame)
m_12 = CSV.read("Data_csv/updated_data/min_12_fix.csv", DataFrame)


fun = CSV.read("Data_csv/updated_data/gumbel_12.csv", DataFrame)

av_0 = CSV.read("Data_csv/updated_data/av_0.csv", DataFrame)
av_12 = CSV.read("Data_csv/updated_data/av_1sqrt2.csv", DataFrame)


########### function to compute the thing

function mu2(mu1, theta2, theta1, lambda::Float64, sigma::Float64)
    results = []
    for k in 1:10
        sum_term = sum(log.(lambda .^(j .-1) .* (theta2 ./ theta1) .^sigma) for j in k:(k+k-1))
        mu2 = mu1 + (1/k) * sum_term
        push!(results, mu2)
    end
    return results
end

###############

g1 = scatter(window_sizes, av_12.scale, xticks=1:1:13,
xlabel = " k ", ylabel = L"\sigma", mc="tomato2",  ms=3, ma=1)
pl.plot!(window_sizes,  av_0.scale[1] ./ window_sizes )

scatter(window_sizes, av_12.location, xticks=1:1:13,
xlabel = " k ", ylabel = L"\mu", mc="tomato2",  ms=3, ma=1)
pl.plot!(window_sizes, av_12.location[1] .+ moving_average.( av_12.location , window_sizes))


updated_mu = mu2.(av_12.location[1], av_12.EI, av_12.EI[1], 2, av_12.scale[1])
av_12.EI
pl.plot!(window_sizes, moving_average.((av_12.location[1]), window_sizes))


log.(2 .^(window_sizes .- 1) .* ( av_12.EI[1] ./ av_12.EI[1] ) .^(av_12.scale[1]))

g2 = scatter(window_sizes, av_12.location, xticks=1:1:13,
xlabel = " k ", ylabel = L"\mu", mc="tomato2",  ms=3, ma=1)
pl.plot!(window_sizes, (av_12.location[1] ./ window_sizes) .+ av_12.location[1] ./ exp.(2))


pl.plot(g1, g2, layout = (1,2), size=(600, 350), title = "independent average")


### Periodic average
g11 = scatter(window_sizes, av_0.scale, xticks=1:1:13,
xlabel = " k ", ylabel = L"\sigma", mc="tomato2",  ms=3, ma=1)
pl.plot!(window_sizes, (av_0.scale[1] ./ window_sizes))

g22 = scatter(window_sizes, av_0.location, xticks=1:1:13,
xlabel = " k ", ylabel = L"\sigma", mc="tomato2",  ms=3, ma=1)
pl.plot!(window_sizes, ( av_0.location[1] ) ./  window_sizes )


### minimums

aaa = scatter(window_sizes, m_0.location, xticks=1:1:13,
xlabel = " k ", ylabel = L"\mu", mc="tomato2",  ms=3, ma=1, legend=false)
pl.plot!(window_sizes, m_0.location[1].+ log.( (m_0.EI ./ m_0.EI[1]).^(m_0.scale[1]))  .- log.(2 .^(window_sizes .- 1)))

savefig(aaa,"Output_Images/updated_plots/periodic_min.pdf")

aaaaaa = scatter(window_sizes, m_0.scale, xticks=1:1:13,
xlabel = " k ", ylabel = L"\sigma", mc="tomato2",  ms=3, ma=1, legend=false)
pl.plot!(window_sizes, m_0.scale[1] ./ 2 .^(window_sizes .-1))


savefig(aaaaaa,"Output_Images/updated_plots/periodic_min_scale.pdf")

function EI_window_min2(orbits, window_sizes)
    location_params = Float64[]
    scale_params = Float64[]
    EI = Float64[]
    shape_params = Float64[]
    
    for windows in window_sizes
        minimums = mapslices(x -> rollmin(x, windows), mat_orb, dims=1)
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

roll_min_matrix = mapslices(x -> rollmin(x, 3), mat_orb, dims=1)


g1 = scatter(window_sizes, testing[2], xticks=1:1:13,
xlabel = " k ", ylabel = L"\sigma", mc="tomato2",  ms=3, ma=1)
pl.plot!(window_sizes,  testing[2][1] ./ window_sizes )

scatter!(window_sizes, testing[1], xticks=1:1:13,
xlabel = " k ", ylabel = L"\mu", mc="tomato2",  ms=3, ma=1)
pl.plot!(window_sizes, (testing[1][1]* 0.9 ./ window_sizes) .+ ( testing[2][1]) ./ testing[4][1])

testing[1]