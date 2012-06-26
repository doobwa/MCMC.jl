function f(x::Array{Float64,2})
  x[1]^2 + x[2]^2
end
function gradient(x::Array{Float64,2})
  [2*x[1] 2*x[2]]
end

# Initial conditions
x0 = [0.0 0.0]

# Tuning parameters
epsilon = .001  # step size
L = 10          # number of steps

# Basic tests on regular interface
x,gx = hmc_sampler(x0, f, gradient, epsilon, L)

@assert isa(x,Array{Float64,2})
@assert isa(gx,Float64)

# Basic tests on Density interface
dd = DifferentiableDensity(f,gradient)
x,gx = hmc_sampler(x0, dd, epsilon, L)

@assert isa(x,Array{Float64,2})
@assert isa(gx,Float64)

# Example of using Options with an HMC sampler allowing bounds
opts = Options(:stepsize,0.01,
               :numsteps,int(50),
               :bounded,[1],          # list of dimensions that are bounded
               :lower_bounds,[0],  # lower bound for 1st dim
               :upper_bounds,[Inf])   # upper bound for 1st dim
x,gx = bounded_hmc_sampler(x0, dd, opts)

@assert isa(x,Array{Float64,2})
@assert isa(gx,Float64)
@assert x[1] > 0
