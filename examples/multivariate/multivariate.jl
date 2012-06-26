# The multivariate density discussed in the following blog posts:
# http://dmbates.blogspot.com/2012_05_01_archive.html
# http://darrenjw.wordpress.com/2011/07/16/gibbs-sampler-in-various-languages-revisited/

function f(v::Array{Float64,2})
  x = v[1]
  y = v[2]
  if x < 0
    return -Inf
  end
  2*log(x) - x*y^2 - y^2 + 2*y - 4*x
end
function gradient(v::Array{Float64,2})
  x = v[1]
  y = v[2]
  # TODO: Deal with x < 0 boundary
  dx = 2/x - y^2 - 4
  dy = -2*x*y - 2*y + 2
  [dx dy]
end

# Sampler that draws from the conditional distributions
function gibbs_sampler(v::Array{Float64}, dd, params)
  x = v[1]
  y = v[2]
  x = rgamma(1,3,1/(y*y + 4))[1]
  y = rnorm(1, 1/(x+1),1/sqrt(2*(x + 1)))[1]
  [x y], dd.f([x y])
end

# General MCMC routine for collecting samples from a chain.  
# We enforce samplers to return the value of the density f at the new state.
function mcmc(x::Array{Float64}, dd, sampler::Function, opts::Options, niter::Int64)
  P = length(x)
  xs = zeros((niter,P))
  progress = zeros(niter)
  for iter = 1:niter
    x, fx = sampler(x, dd, opts)
    xs[iter,:] = x
    progress[iter] = fx
  end
  return xs, progress
end

x0 = [0.1 0.5]
num_iterations = 10000

dd = DifferentiableDensity(f,gradient)

opts = Options(:stepsize,0.0001,
               :numsteps,int(50),
               :bounded,[1],          # list of dimensions that are bounded
               :lower_bounds,[0],     # lower bound for 1st dim
               :upper_bounds,[Inf])   # upper bound for 1st dim
xs, progress = mcmc(x0, dd, bounded_hmc_sampler, opts, num_iterations)
csvwrite("examples/multivariate/hmc.dat",xs)

opts = Options()
xs, progress = mcmc(x0, dd, gibbs_sampler, opts, num_iterations)
csvwrite("examples/multivariate/gibbs.dat",xs)

opts = Options()
xs, progress = mcmc(x0, dd, mv_slice_sampler, opts, num_iterations)
csvwrite("examples/multivariate/slice.dat",xs)
