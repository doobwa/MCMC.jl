module MCMC
  using Base

  export hmc_sampler, bounded_hmc_sampler, mh_sampler, slice_sampler

  #load("Rmath.jl")

  load("Distributions")
  using Distributions

  load("Options")
  using OptionsMod

  load("MCMC/src/density.jl")
  load("MCMC/src/mhsampler.jl")
  load("MCMC/src/slicesampler.jl")
  load("MCMC/src/hmcsampler.jl")
end
