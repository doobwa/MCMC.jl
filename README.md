# MCMC Routines in Julia

In Bayesian modeling we are often faced with an unnormalized density (e.g. a posterior distribution) that we want to know about.  A few general purpose MCMC techniques for such situations include [Metropolis-Hastings](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm), [slice sampling](https://en.wikipedia.org/wiki/Slice_sampling), and [Hamiltonian Monte Carlo](http://www.cs.toronto.edu/~radford/ftp/ham-mcmc.pdf) (HMC).  

By implementing these samplers in [Julia](http://julialang.org), they can be written naturally while being quite fast.  By providing a common interface across samplers, we can perform MCMC with much less effort: in several of the examples below we only needed to implement a function that computes the unnormalized density of interest.  More often, these tools are  incorporated as components of a larger sampling scheme.

## Basic idea

Suppose we have a function `g(x)`.  Our aim is to treat `g` as an unnormalized probability density and obtain samples from the resulting probability distribution.  We can use slice sampling to obtain a new sample from `g` where the previous state was `x` using the `slice_sampler` method in this package:

    x,gx = slice_sampler(x,g)

This function also returns `g(x)` at the new value, as shown above.  We can obtain many samples of by repeatedly calling the above function.  In truth, this sequence of `x` values is a Markov chain whose limiting distribution is `g(x)`, and though the samples are not independent, they often can be a useful summary of the function `g(x)`.  In the case of Bayesian modeling, the function of interest `g(x)` is the posterior distribution of a model's parameters given the observed data.

## Univariate example

Consider the following two functions, where `dnorm(x,mu,sigma)` is the density Normal(x; mu,sigma):

    g(x) = dnorm(x,0,.5)
    h(x) = dnorm(x,-1,.35) + 2 * dnorm(x,1,.5) + 1.5 * dnorm(x,2,.25)

The following will run a julia script that performs 10000 iterations of slice sampling and MH on `log(g)` and `log(h)` respectively.  From the `mcmc/` directory:

    julia> load("examples/example.jl")

On my laptop I get the following timings:
    slice on g: 30.5 ms
    slice on h: 90.0 ms
    mh on g: 9.2 ms
    mh on h: 15.3 ms

An R implementation of slice sampling is included for comparison (in seconds):

    > source("examples/example.r")
    > system.time(samples <- mcmc(0,lg,slice,niter=10000))
       user  system elapsed 
      1.420   0.000   1.421 
    > system.time(samples <- mcmc(0,lh,slice,niter=10000))
       user  system elapsed 
      2.050   0.000   2.041 

For slice sampling `h` we see about a 18x speedup.  This is encouraging because in many applications it is expensive to repeatedly compute the density, but with these functions we can now write them in Julia.

Both the Julia version and the R version are using R libraries for random uniform and exponential draws during the slice sampling procedure.

The R script `examples/plot.r` loads the various samples and plots summaries.  [This image](https://github.com/doobwa/mcmc.jl/tree/master/examples/results/slice.h.dat.png) shows the results of using the slice sampler in Julia on `log(h(x))`.  The true density is the black curve superimposed on the histogram.  These plots provide reassurance that we are correctly sampling from the intended distribution and mixing fairly well.

Note: The goal here is to illustrate the use of these routines, not compare the efficacy of different samplers.  

## Multivariate example

In [several](http://darrenjw.wordpress.com/2010/04/28/mcmc-programming-in-r-python-java-and-c/) [recent](https://darrenjw.wordpress.com/2011/07/31/faster-gibbs-sampling-mcmc-from-within-r/) [blog](http://dirk.eddelbuettel.com/blog/2011/07/14/) [posts](http://dmbates.blogspot.com/2012_05_01_archive.html) people have explored MCMC for the following bivariate density:

f(x,y) = K x^2 exp{-xy^2 - y^2 + 2y - 4x}

where K is some unknown normalizing constant.  These posts perform MCMC by [Gibbs sampling](https://en.wikipedia.org/wiki/Gibbs_sampling) using the following full conditional distributions:

x|y ~ Gamma(3,1/(y^2+4))
y|x ~ Normal(1/(1+x), 1/(2(1+x)))

[NB: The second parameter in the Gamma is a rate parameter.]
For some models, it can be a pain (or analytically intractable) to derive full conditional distributions for some parameters.  The `examples/multivariate/multivariate.jl` script compares the above approach versus iteratively slice sampling in each dimension.  [These plots](https://github.com/doobwa/mcmc.jl/tree/master/examples/multivariate/compare.png) show the density evaluated at a grid of points, then 2000 samples from the Gibbs sampling algorithm, then 2000 samples from the multivariate slice sampler.  One can also look at the [trace plots](https://github.com/doobwa/mcmc.jl/tree/master/examples/multivariate/trace.png) for each sampler across each iteration.

Here is an excerpt from the example, showing how we are calling the `mcmc` method with a particular `DifferentiableDensity` object `dd`.  When using the `gibbs_sampler` we can obtain 10000 samples in .043 seconds.  Slice sampling here took about 3 times longer, but remember this is quite encouraging as this method is using the default parameters, it's simply sampling each dimension in turn, and required no special knowledge about the density `f(x,y)`. 

```
julia> @time xs, progress = mcmc(x0, dd, gibbs_sampler, opts, num_iterations);
elapsed time: 0.04285693168640137 seconds

julia> @time xs, progress = mcmc(x0, dd, mv_slice_sampler, opts, num_iterations);
elapsed time: 0.12020492553710938 seconds
```

## Credits
The slice sampling and HMC algorithms implemented here are based on Radford Neal's [R code](http://www.cs.toronto.edu/~radford/software-online.html).

## TODO
* Implement Radford Neal's...
  - windowed HMC
  - tempered HMC
  - [x] bounded HMC
* Debug HMC
* Allow samplers to return other metadata (e.g. number of function evaluations)
* Consistently use the previous value of `g(x)` so that we do fewer evaluations of `g`
* [x] Make use of Density type
* [x] Use Options type for passing options to the various samplers
* Add more checks for sensible parameter settings in each sampler method
