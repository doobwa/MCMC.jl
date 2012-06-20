function mh_sampler(x::Float64, g::Function, sd::Float64, gx::Float64)
  x1 = x + randn()*sd
  gx1 = g(x1)
  if gx1 - gx > log(rand())
    x = x1
    gx = gx1
  end
  return x,gx
end

function mh_sampler(x::Float64, g::Function, sd::Float64)
  gx = g(x)
  mh_sampler(x,g,sd,gx)
end

mh_sampler(x::Float64, g::Function) = mh_sampler(x,g,.5)

# Interface using Density type
mh_sampler(x::Float64, g::Density, sd::Float64, gx::Float64) = mh_sampler(x, Density.f, sd, gx)

