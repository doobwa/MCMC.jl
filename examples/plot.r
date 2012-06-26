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

# Plots for multivariate example
samples_gibbs <- read.csv("examples/multivariate/gibbs.dat",header=FALSE)
samples_hmc <- read.csv("examples/multivariate/hmc.dat",header=FALSE)
samples_slice <- read.csv("examples/multivariate/slice.dat",header=FALSE)
png("compare.png",width=1000,height=300)
par(mfrow=c(1,3))
plot(samples_gibbs)
plot(samples_hmc)
plot(samples_slice)
  acf(samples)
  dev.off()
}
