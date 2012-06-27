load("src/mcmc.jl")

function g(x)
  log(dnorm(x,0,.5))
end

function h(x)
  log(dnorm(x,-1,.35) + 2*dnorm(x,1,.5) + 1.5*dnorm(x,2,.25))
end

function mcmc(x::Float64, g::Function, sampler::Function, niter::Int64)
  xs = zeros(niter)
  gx = g(x)
  for iter = 1:niter
    x,gx = sampler(x,g,gx)
    xs[iter] = x
  end
  return xs
end

x0 = 0.0
niter = 10000

@elapsed xs = mcmc(x0,g,slice_sampler,niter)
csvwrite("examples/results/slice.g.dat",xs)

@elapsed xs = mcmc(x0,h,slice_sampler,niter) 
csvwrite("examples/results/slice.h.dat",xs)

@elapsed xs = mcmc(x0,g,mh_sampler,niter)
csvwrite("examples/results/mh.g.dat",xs)

@elapsed xs = mcmc(x0,h,mh_sampler,niter)
csvwrite("examples/results/mh.h.dat",xs)

