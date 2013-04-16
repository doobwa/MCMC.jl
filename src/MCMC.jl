require("Distributions")
require("Options")
  
module MCMC

  using Base
  using Distributions
  using OptionsMod

  export hmc_sampler, bounded_hmc_sampler, mh_sampler, slice_sampler

  #load("Rmath.jl")

  include("density.jl")
  include("mhsampler.jl")
  include("slicesampler.jl")
  include("hmcsampler.jl")
end
