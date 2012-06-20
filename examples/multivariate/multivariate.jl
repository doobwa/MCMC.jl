# The multivariate density discussed in the following blog posts:
# http://dmbates.blogspot.com/2012_05_01_archive.html
# http://darrenjw.wordpress.com/2011/07/16/gibbs-sampler-in-various-languages-revisited/

function f(v::Array{Float64,2})
  x = v[1]
  y = v[2]
  x^2 * exp(-x*y^2 - y^2 + 2*y - 4*x)
end
function gradient(v::Array{Float64,2})
  x = v[1]
  y = v[2]
  tmp = exp(-x*y^2 - y^2 + 2*y - 4*x)
  dx = 2 * x * tmp + x^2 * tmp * (-2*x*y - 4)
  dy = x^2 * tmp * (-2*x*y - 2*y + 2)
  [dx dy]
end

f([.2 .5])
gradient([.2 .5])
dd = DifferentiableDensity(f,gradient)

# Sampler that draws from the conditional distributions
function gibbs_sampler(v::Array{Float64}, dd, params)
  x = v[1]
  y = v[2]
  x = randg(3) * (y*y + 4)
  y = 1/(x + 1) + randn()/sqrt(2(x + 1))
  [x y], dd.f([x y])
end

# General MCMC routine for collecting samples from a chain.  
# We enforce samplers to return the value of the density f at the new state.
function mcmc(x::Array{Float64}, dd::DifferentiableDensity, sampler::Function, sampler_params::Array, niter::Int64)
  P = length(x)
  xs = zeros((niter,P))
  progress = zeros(niter)
  for iter = 1:niter
    x, gx = sampler(x, dd, sampler_params)
    xs[iter,:] = x
    progress[iter] = gx
  end
  return xs, progress
end

x0 = [0.0 0.0]
hmc_params = [0.05 20]
num_iterations = 10000
xs, progress = mcmc(x0, dd, hmc_sampler, hmc_params, num_iterations)
csvwrite("examples/multivariate/hmc.dat",xs)

gibbs_params = [3 4 1]  # TODO: Use this in the actual method
xs, progress = mcmc(x0, dd, gibbs_sampler, gibbs_params, num_iterations)
csvwrite("examples/multivariate/gibbs.dat",xs)


# Function for Gibbs sampling multivariate density via slice sampling
# TODO: Following is still broken
function mv_slice_sampler(x::Array{Float64}, f::Function)
  P = length(x)
  function fp(w)  # want to fix all but the p'th value of this array
    x[p] = w
    f(x)
  end
  for p = 1:length(x)
    x[p], gxp = slice_sampler(x[p], fp)  # TODO: initialize at last point
  end
  x, gxp
end
