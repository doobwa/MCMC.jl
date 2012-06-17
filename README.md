# MCMC Routines in Julia

In Bayesian modeling we are often faced with an unnormalized density (e.g. a posterior distribution) that we want to know about.  A few general purpose MCMC techniques for such situations include [Metropolis-Hastings](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm), [slice sampling](https://en.wikipedia.org/wiki/Slice_sampling), and [Hamiltonian Monte Carlo](http://www.cs.toronto.edu/~radford/ftp/ham-mcmc.pdf) (HMC).

By implementing these samplers in (Julia)[http://julialang.org], they can be written naturally while being quite fast.  By providing a common interface across samplers, it also becomes easy to incorporate these techniques as components of larger projects.

## Example

Use slice sampling to obtain a new sample `x` from `g` where the previous state was `x`.  Also returns `g(x)` for the new value.

    x,gx = slice_sampler(x,g)

## Comparison with R

Consider the following two functions, where `dnorm(x,mu,sigma)` is the density Normal(x; mu,sigma):

    g(x) = dnorm(x,0,.5)
    h(x) = dnorm(x,-1,.35) + 2 * dnorm(x,1,.5) + 1.5 * dnorm(x,2,.25)

The following will run a julia script that performs 10000 iterations of slice sampling and MH on log(g) and log(h) respectively.  From the `mcmc/` directory:

    julia> load("[path to julia]/extras/Rmath.jl")
    julia> load("examples/example.jl")
    slice on g took 78.9 ms
    slice on h took 146.0 ms
    mh on g took 15.8 ms
    mh on h took 12.3 ms

An R implementation of slice sampling is included for comparison.

    > source("examples/example.r")
    > system.time(samples <- mcmc(0,lg,slice,niter=10000))
       user  system elapsed 
      1.420   0.000   1.421 
    > system.time(samples <- mcmc(0,lh,slice,niter=10000))
       user  system elapsed 
      2.050   0.000   2.041 

The R script `examples/plot.r` loads the various samples and plots summaries.  (This image)[https://github.com/doobwa/mcmc.jl/tree/master/examples/results/slice.h.dat.png] shows the results of using the slice sampler in Julia on `log(h(x))`.  The true density is the black curve superimposed on the histogram.  These plots provide reassurance that we are correctly sampling from the intended distribution and mixing fairly well.

Note: The goal here is to illustrate the use of these routines, not compare the efficacy of different samplers.  

The slice sampling code is due to Radford Neal, obtained (here)[http://www.cs.toronto.edu/~radford/software-online.html].
