load("Distributions")
using Distributions
load("MCMC")
using MCMC

load("benchmarks/julia/densities.jl")

function mcmc(log_density::Function, x::Vector, sampler::Function, n_chains::Int64, n_samples::Int64)
    xs = zeros(n_samples, length(x))
    gx = log_density(x)
    for iter = 1:n_samples
        x, gx = sampler(log_density, x, gx)
        xs[iter,:] = x
    end
    return xs
end

function multivariate_slice(f::Function, x::Vector, sampler::Function)
    P = length(x)
    fx = f(x)
    for p = 1:length(x)
        # Keep all but the p'th value fixed
        function fp(w)  
            x[p] = w
            f(x)
        end
        
        # Pass this to a univariate slice sampler
        x[p], fx = slice_sampler(x[p], fp, fx)
    end
    x, fx
end
mv_slice(f::Function, x, fx) = multivariate_slice(f, x, slice_sampler)
mv_slice(f::Function, x, fx) = multivariate_slice(f, x, slice_sampler)

function run_australian()
    d = csvread("benchmarks/csv/australian.csv")
    y = d[:,1]
    x = d[:,2:end]
    theta = rand(MultivariateNormal([0. for i in 1:14], .000001 *  eye(14)))
    theta[13] = 0
    samples = mcmc(logistic_regression, theta, mv_slice, 1, 10)
end
