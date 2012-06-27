source("examples/example.r")

# Compare true density to samples collected from an MCMC chain.
# Also show the samples by iteration and the autocorrelation function.
#
# g: true density function
# samples: vector of samples from a density g
plot.chain <- function(g,samples) {
  par(mfrow=c(1,3))#,mar=c(2,1,1,1))
  window <- .01
  xs <- seq(min(samples),max(samples),window)
  hist(samples, breaks=20, probability=TRUE,main="")
  lines(xs,g(xs)/sum(g(xs))/window,lwd=2)
  plot(samples, type="l")
  acf(samples)
}

# Plot of the samples from the R code and the julia code.
fns <- list(g,h)
for (i in c("g","h")) {   # target densities
  chains <- c(paste("examples/results/r.slice.",i,".dat",sep=""),
              paste("examples/results/slice.",i,".dat",sep=""),
              paste("examples/results/mh.",i,".dat",sep=""))
  f <- get(i)
  for (chain in chains) { # methods
    samples <- read.csv(chain)[,1]
    png(paste(chain,".png",sep=""),width=1000,height=300)
    plot.chain(f,samples)
    dev.off()
  }
}

# Plots for multivariate example.
f <- function(x,y) x^2 * exp(-x*y^2 - y^2 + 2*y - 4*x)
xs <- seq(0,3,by=.01)
ys <- seq(-1,2,by=.01)
truth <- matrix(0,length(xs),length(ys))
for (i in 1:length(xs)) {
  for (j in 1:length(ys)) {
    truth[i,j] <- f(xs[i],ys[j])
  }
}

samples_gibbs <- read.csv("examples/multivariate/gibbs.dat",header=FALSE)
samples_hmc <- read.csv("examples/multivariate/hmc.dat",header=FALSE)
samples_slice <- read.csv("examples/multivariate/slice.dat",header=FALSE)
png("examples/multivariate/compare.png",width=1000,height=300)
xlims <- c(0,3)
ylims <- c(-1,2)
par(mfrow=c(1,3))
image(xs,ys,truth)
plot(samples_gibbs[1:1000,],xlim=xlims,ylim=ylims,xlab="Gibbs",ylab="")
plot(samples_slice[1:1000,],xlim=xlims,ylim=ylims,xlab="Slice",ylab="")
#plot(samples_hmc[1:1000,],xlim=xlims,ylim=ylims,xlab="HMC",ylab="")
#acf(samples_slice)
dev.off()

png("examples/multivariate/trace.png",width=500,height=500)
par(mfrow=c(2,2))
plot(samples_gibbs[1:1000,1],type="l",xlab="Iteration",ylab="Gibbs sampling x",ylim=xlims)
plot(samples_slice[1:1000,1],type="l",xlab="Iteration",ylab="Slice sampling x",ylim=xlims)
plot(samples_gibbs[1:1000,2],type="l",xlab="Iteration",ylab="Gibbs sampling y",ylim=ylims)
plot(samples_slice[1:1000,2],type="l",xlab="Iteration",ylab="Slice sampling y",ylim=ylims)
dev.off()
