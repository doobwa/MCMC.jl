# Function for Gibbs sampling multivariate density via slice sampling
function mv_slice_sampler(x::Array{Float64}, f::Function, opts::Options)
  P = length(x)
  gxp = 0
  for p = 1:length(x)
    function fp(w)  # want to fix all but the p'th value of this array
      x[p] = w
      f(x)
    end
    # TODO: initialize at last point
    # TODO: Pass on options to univariate slice sampler
    x[p], gxp = slice_sampler(x[p], fp)  
  end
  x, gxp
end
mv_slice_sampler(x::Array{Float64}, d::Union(Density,DifferentiableDensity), opts::Options) = mv_slice_sampler(x, d.f, opts)

