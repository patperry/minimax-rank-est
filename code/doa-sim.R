# doa-sim.R
# ---------
source("doa-sim-params.R")

SampleDOASim <- function(par.sim) {
  # Sample a DOA simulation object
  #
  # Args:
  #   par.sim: a DOASimParams object with the simulation parameters
  #
  # Returns:
  #   a `DOASim' object sampled from the given simulation
  if (!inherits(par.sim, "DOASimParams")) {
    stop("Argument `par.sim' must be a `DOASimParams' object, not: ",
         class(par.sim), ".") 
  }
  
  num.snapshots <- par.sim$num.times
  dim           <- par.sim$dim

  signal   <- SampleCovTimeSeries(par.sim$cov.signal)
  noise    <- matrix(complex(real=rnorm(num.snapshots * dim, sd=sqrt(1/2)),
                             imag=rnorm(num.snapshots * dim, sd=sqrt(1/2))),
                     num.snapshots, dim)
  snapshot <- signal + noise

  res <- list(par=par.sim,
    num.snapshots=num.snapshots,
              dim=dim,
             time=par.sim$time,
         snapshot=snapshot,
           signal=signal,
            noise=noise) 
  class(res) <- c("DOASim", class(res))
  res
}

