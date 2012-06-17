source("examples/slice.r")

mcmc <- function(x,g,method,niter=1000) {
  samples <- rep(0,niter)
  for (iter in 1:niter) {
    x <- method(x,g)
    samples[iter] <- x
  }
  return(samples)
}

g <- function(x) dnorm(x,0,.5)
h <- function(x) dnorm(x,-1,.35) + 2*dnorm(x,1,.5) + 1.5*dnorm(x,2,.25)

lg <- function(x) log(g(x))
lh <- function(x) log(h(x))

slice <- function(x,g) {
  uni.slice.alt(x,g,w=1,m=1000)
}

# Time 10000 iterations from slice sampler and save results for g and h
print(system.time(samples <- mcmc(0,lg,slice,niter=10000)))
write.csv(samples,file="examples/results/r.slice.g.dat",row.names=FALSE)

print(system.time(samples <- mcmc(0,lh,slice,niter=10000)))
write.csv(samples,file="examples/results/r.slice.h.dat",row.names=FALSE)

