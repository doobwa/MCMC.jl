# Adapted directly from Radford Neal's R code

# Alternative interface using Density type
hmc_sampler(x::Array{Float64}, dd::DifferentiableDensity, epsilon::Float64, L::Int64) = hmc_sampler(x, dd.f, dd.gradient, epsilon, L)
# Some default values for epsilon and L
hmc_sampler(x::Array{Float64}, dd::DifferentiableDensity, gx::Float64) = hmc_sampler(x,dd,.01,20)

function hmc_sampler(x::Array{Float64}, dd::DifferentiableDensity, params::Array) 
  L = convert(Int64, params[2])
  hmc_sampler(x,dd,params[1],L)  # TODO Use a different data structure for params
end

function hmc_sampler(current_q::Array{Float64}, U::Function, grad_U::Function, epsilon::Float64, L::Int64)

  q = current_q
  p = randn(length(q))'  # independent standard normal variates (row vector)
  current_p = p

  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2

  # Alternate full steps for position and momentum
  for i in 1:L
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if i!=L
        p = p - epsilon * grad_U(q)
    end
  end

  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2

  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p

  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p.^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p.^2) / 2

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if rand() < exp(current_U-proposed_U+current_K-proposed_K)
    return q, proposed_U  # accept
  else
    return current_q, current_U  # reject
  end
end
