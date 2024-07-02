using Distributions
using Plots

σ = 2/3
μ = 1/2

θ_values = 0.1:0.1:1

#μ_k = μ + σ*log(θ)
# Generate a range of values for the x-axis
x = -5:0.01:5

# Create the animation
anim = @animate for θ in θ_values
    d = Gumbel(μ + σ*log(θ), σ)
    y = pdf.(d, x)
    plot(x, y, title="Gumbel Distribution", xlabel="x", ylabel="pdf(x)", legend=false, ylim=(0, 1))
end

# Save the animation as a gif
gif(anim, "gumbel_distribution.gif", fps=10)


n = 1/2
k = 1.2

function aaaaa(n, k)
    test = (log(n))^(-n/k) == log(k/ ( (-(k/n)*log(n))^(1+n/k) ))
    return test
end

aaaaa(n,k)