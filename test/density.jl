function f(x::Array{Float64,2})
  x[1]^2 + x[2]^2
end
function gradient(x::Array{Float64,2})
  [2*x[1] 2*x[2]]
end

x = [1.0 3.0]
dd = DifferentiableDensity(f,gradient)
@assert dd.f(x) == 10.0
@assert dd.gradient(x) == [2.0 6.0]
