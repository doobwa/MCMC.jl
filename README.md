# MCMC Routines in Julia

In Bayesian modeling we are often faced with an unnormalized density (e.g. a posterior distribution) that we want to know about.  A few general purpose MCMC techniques for such situations include [Metropolis-Hastings](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm), [slice sampling](https://en.wikipedia.org/wiki/Slice_sampling), and [Hamiltonian Monte Carlo](http://www.cs.toronto.edu/~radford/ftp/ham-mcmc.pdf) (HMC).  

By implementing these samplers in [Julia](http://julialang.org), they can be written naturally while being quite fast.  By providing a common interface across samplers, we can perform MCMC with much less effort: in several of the examples below we only needed to implement a function that computes the unnormalized density of interest.  More often, these tools are  incorporated as components of a larger sampling scheme.

## Example

Suppose we have a function `g(x)` that we want to sample.  We can use slice sampling to obtain a new sample `x_new` from `g` where the previous state was `x`.  The sampler also returns `g(x)` for the new value.

    x_new,gx = slice_sampler(x,g)

## Comparison with R

Consider the following two functions, where `dnorm(x,mu,sigma)` is the density Normal(x; mu,sigma):

    g(x) = dnorm(x,0,.5)
    h(x) = dnorm(x,-1,.35) + 2 * dnorm(x,1,.5) + 1.5 * dnorm(x,2,.25)

The following will run a julia script that performs 10000 iterations of slice sampling and MH on `log(g)` and `log(h)` respectively.  From the `mcmc/` directory:

    julia> load("[path to julia]/extras/Rmath.jl")
    julia> load("examples/example.jl")
    slice on g took 78.9 ms
    slice on h took 146.0 ms
    mh on g took 15.8 ms
    mh on h took 12.3 ms

An R implementation of slice sampling is included for comparison (in seconds):

    > source("examples/example.r")
    > system.time(samples <- mcmc(0,lg,slice,niter=10000))
       user  system elapsed 
      1.420   0.000   1.421 
    > system.time(samples <- mcmc(0,lh,slice,niter=10000))
       user  system elapsed 
      2.050   0.000   2.041 

For slice sampling `h` we see about a 14x speedup.  This is encouraging because in many applications it is expensive to repeatedly compute the density, but with these functions we can now write them in Julia.

Both the Julia version and the R version are using R libraries for random uniform and exponential draws during the slice sampling procedure.

The R script `examples/plot.r` loads the various samples and plots summaries.  [This image](https://github.com/doobwa/mcmc.jl/tree/master/examples/results/slice.h.dat.png) shows the results of using the slice sampler in Julia on `log(h(x))`.  The true density is the black curve superimposed on the histogram.  These plots provide reassurance that we are correctly sampling from the intended distribution and mixing fairly well.

Note: The goal here is to illustrate the use of these routines, not compare the efficacy of different samplers.  

## Credits
The slice sampling code is due to Radford Neal, obtained [here](http://www.cs.toronto.edu/~radford/software-online.html).

## TODO
* Before adding in HMC, it might be nice to use Julia's Type system to allow for the presence or absence of a gradient for the function of interest.
* Use the previous value of `g(x)` so that we do fewer evaluations of `g`
