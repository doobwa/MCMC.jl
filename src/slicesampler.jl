# Ported from Radford Neal's R code, with a few thinggs missing

# Arguments:
#
#   x0    Initial point
#   g     Function returning the log of the probability density (plus constant)
#   w     Size of the steps for creating interval (default 1)
#   m     Limit on steps
#   lower Lower bound on support of the distribution (default -Inf)
#   upper Upper bound on support of the distribution (default +Inf)
#   gx0   g(x0)
#
function slice_sampler(x0::Float64, g::Function, w::Float64, m::Int64, lower::Float64, upper::Float64, gx0::Float64)

  if w <= 0
    error("Negative w not allowed")
  end
  if m <= 0
    error("Limit on steps must be positive")
  end
  if upper < lower
    error("Upper limit must be above lower limit")
  end
  
  # Determine the slice level, in log terms.
  
  logy = gx0 - rand(Exponential(1.0))

  # Find the initial interval to sample from.

  u = rand() * w
  L = x0 - u
  R = x0 + (w-u)

  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.  

  J = floor(rand() * m)
  K = (m-1) - J

  while J > 0
    if L <= lower || g(L) <= logy
      break
    end
    L -= w
    J -= 1
  end

  while K > 0
    if R >= upper || g(R) <= logy
      break
    end
    R += w
    K -= 1
  end

  # Shrink interval to lower and upper bounds.

  L = L < lower ? lower : L
  R = R > upper ? upper : R
  
  # Sample from the interval, shrinking it on each rejection.

  x1 = 0  # need to initialize it in this scope first
  gx1 = 0
  while true 
    x1 = rand() * (R-L) + L
    gx1 = g(x1)
    if gx1 >= logy
      break
    end
    if x1 > x0
      R = x1
    else
      L = x1
    end
  end

  return x1,gx1
end

function slice_sampler(x0::Float64, g::Function, w::Float64, m::Int64, lower::Float64, upper::Float64)
  gx0 = g(x0)
  slice_sampler(x0,g,w,m,lower,upper,gx0)
end

function slice_sampler(x0::Float64, g::Function, gx0::Float64)
  slice_sampler(x0,g,.5,10000,-Inf,Inf,gx0)
end

function slice_sampler(x0::Float64, g::Function)
  slice_sampler(x0,g,.5,10000,-Inf,Inf)
end

# Interface using Density type

slice_sampler(x0::Float64, d::Density, w::Float64, m::Int64, lower::Float64, upper::Float64, gx0::Float64) = slice_sampler(x0::Float64, d.f, w::Float64, m::Int64, lower::Float64, upper::Float64, gx0::Float64)

slice_sampler(x0::Float64, d::Density, gx0::Float64) = slice_sampler(x0::Float64, d.f, gx0::Float64)

function slice_sampler(x0::Float64, d::Density, opts::Options)
  gx0 = g(x0)
  slice_sampler(x0::Float64, d.f, opts[:w],opts[:m],opts[:lower],opts[:upper],gx0::Float64)
end
