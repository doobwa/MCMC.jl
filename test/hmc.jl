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
x,gx = hmc_sampler(x0, dd, epsilon, L)

@assert isa(x,Array{Float64,2})
@assert isa(gx,Float64)

