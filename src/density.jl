# Density type that is useful across a variety of the sampling methods.
# NB: The density computed by f need not be normalized.

type Density
  f::Function
  Density(f::Function) = new(f)
end

type DifferentiableDensity  # maybe <: DifferentiableFunction ?
  f::Function 
  gradient::Function
  DifferentiableDensity(f::Function,gradient::Function) = new(f,gradient)
end
