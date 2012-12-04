
## Univariate case
function metropolis_sampler(x::Float64, g::Function, sd::Float64, gx::Float64)
  x1 = x + randn()*sd
  gx1 = g(x1)
  if gx1 - gx > log(rand())
    x = x1
    gx = gx1
  end
  return x,gx
end

function metropolis_sampler(x::Float64, g::Function, sd::Float64)
  gx = g(x)
  metropolis_sampler(x,g,sd,gx)
end

metropolis_sampler(x::Float64, g::Function) = metropolis_sampler(x,g,.5)

# Interface using Density type
metropolis_sampler(x::Float64, g::Density, sd::Float64, gx::Float64) = metropolis_sampler(x, Density.f, sd, gx)

## Multivariate case
function metropolis_sampler(x::Vector{Float64}, g::Function, sd::Float64, gx::Vector{Float64})
  x1 = x + sd*randn(numel(x))
  gx1 = g(x1)
  if gx1 - gx > log(rand())
    x = x1
    gx = gx1
  end
  return x,gx
end

function metropolis_sampler(x::Vector{Float64}, g::Function, sd::Matrix{Float64}, gx::Vector{Float64})
# If C is the desired covariance of the proposal then sd = chol(C)'
# The reason for pre-multiplying is to generate a vector rather than a 
# matrix.
  x1 = x + sd*randn(numel(x))
  gx1 = g(x1)
  if gx1 - gx > log(rand())
    x = x1
    gx = gx1
  end
  return x,gx
end

# Add in versions with default arguments
function metropolis_sampler(x::Vector{Float64}, g::Function, sd::Matrix{Float64})
  gx = g(x)
  metropolis_sampler(x,g,sd,gx)
end

metropolis_sampler(x::Float64, g::Function) = metropolis_sampler(x,g,.5*eye(numel(x)))

# Interface using Density type
metropolis_sampler(x::Vector{Float64}, g::Density, sd::Matrix{Float64}, gx::Vector{Float64}) = metropolis_sampler(x, Density.f, sd, gx)
