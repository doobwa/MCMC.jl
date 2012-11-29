# Univariate case
function metrop_sampler(x::Float64, g::Function, sd::Float64, gx::Float64)
  x1 = x + randn()*sd
  gx1 = g(x1)
  if gx1 - gx > log(rand())
    x = x1
    gx = gx1
  end
  return x,gx
end

function metrop_sampler(x::Float64, g::Function, sd::Float64)
  gx = g(x)
  mh_sampler(x,g,sd,gx)
end

metrop_sampler(x::Float64, g::Function) = mh_sampler(x,g,.5)

# Interface using Density type
metrop_sampler(x::Float64, g::Density, sd::Float64, gx::Float64) = mh_sampler(x, Density.f, sd, gx)

# Multivariate case
function metrop_sampler(x::Array{Float64,1}, g::Function, sd::Float64, gx::Array{Float64,1})
  x1 = x + sd*randn(numel(x))
  gx1 = g(x1)
  if gx1 - gx > log(rand())
    x = x1
    gx = gx1
  end
  return x,gx
end

function metrop_sampler(x::Vector{Float64}, g::Function, sd::Matrix{Float64}, gx::Vector{Float64})
# If C is the desired covariance of the proposal then sd = chol(C)
  x1 = x + sd*randn(1,numel(x))
  gx1 = g(x1)
  if gx1 - gx > log(rand())
    x = x1
    gx = gx1
  end
  return x,gx
end

# Add in versions with default arguments
