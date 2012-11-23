load("Distributions")
using Distributions

load("MCMC")
using MCMC

function f(x)
  log(0.75 * pdf(Normal(0.0, 1.0), x) + 0.25 * pdf(Normal(10.0, 1.0), x))
end

function g(x)
  log(pdf(Normal(x,0,.5)))
end

function h(x)
  log(pdf(Normal(x,-1,.35)) + 2*pdf(Normal(x,1,.5)) + 1.5*pdf(Normal(x,2,.25)))
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

@elapsed samples = mcmc(0.0, f, slice_sampler, 50_000)
csvwrite("examples/results/slice.f.csv", samples)

@elapsed samples = mcmc(x0,g,slice_sampler,niter)
csvwrite("examples/results/slice.g.csv",samples)

@elapsed samples = mcmc(x0,h,slice_sampler,niter) 
csvwrite("examples/results/slice.h.csv",samples)

@elapsed samples = mcmc(x0,g,mh_sampler,niter)
csvwrite("examples/results/mh.g.csv",samples)

@elapsed samples = mcmc(x0,h,mh_sampler,niter)
csvwrite("examples/results/mh.h.csv",samples)

