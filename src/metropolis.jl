
# Univariate case
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

# Multivariate case
function metropolis_sampler(x::Vector{Float64}, g::Function, sd::Float64, gx::Vector{Float64})
  x1 = x + randn(numel(x))*sd
  gx1 = g(x1)
  if gx1 - gx > log(rand())
    x = x1
    gx = gx1
  end
  return x,gx
end

function metropolis_sampler(x::Vector{Float64}, g::Function, C::Matrix{Float64}, gx::Vector{Float64})
# If C is the desired covariance of the proposal then sd = chol(C)'
# Should probably add some error checking inc ase C is not p.s.d.
  #sd = chol(C)' # Probably should pass sd in directly
  #x1 = x + sd*randn(numel(x))
  x1 = rand(MultivariateNormal(x,C))
  gx1 = g(x1)
  if gx1 - gx > log(rand())
    x = x1
    gx = gx1
  end
  return x,gx
end

# Add in versions with default arguments
function metropolis_sampler(x::Vector{Float64}, g::Function, C::Matrix{Float64})
  gx = g(x)
  metropolis_sampler(x,g,C,gx)
end

metropolis_sampler(x::Float64, g::Function) = metropolis_sampler(x,g,.5*eye(numel(x)))

# Interface using Density type
metropolis_sampler(x::Vector{Float64}, g::Density, C::Matrix{Float64}, gx::Vector{Float64}) = metropolis_sampler(x, Density.f, C, gx)
