# rank-est-generic.R
# ------------------

source("cov-time-series.R")

EstimateRankDummy <- function(num.snapshots, dim, num.signals.0, noise.sd.0,
                              evalue, evalue.sqrt, evector, ...) {
  list(num.signals=num.signals.0, noise.sd=noise.sd.0)
}

MakeRankEstTimeSeriesFun <- function(estimate) {
  # Take a rank estimation rule and return a function that applies the
  # rule to an estimated covariance time series.
  #
  # Args:
  #   estimate: a function that takes arguments `num.snapshots', `dim',
  #             `num.signals.0', `noise.sd.0', `evalues', `evalues.sqrt',
  #             and `evectors'; the function should return a list with
  #             components `num.signals' and `noise.sd'; see
  #             `EsttimateRankDummy'
  #
  # Returns:
  #   a function that applies the rule to a covariance estimate time series  
  function(cov.est.ts, num.signals.0=0, noise.sd.0=1, ...) {
    # Take an estimated covariance time series and estimate the low-rank
    # signal part.
    #
    # Args:
    #   cov.est.ts:    a `CovEstTimeSeries' object with the estimated time-
    #                  series
    #   num.signals.0: the initial estimate for the number of signals
    #   noise.sd.0:    the initial estimate for the standard deviation of the
    #                  noise
    #
    # Returns:
    #   A list with components
    #     num.times: the number of time points in the time-series
    #     time:      the times
    #     cov.ts:    an `CovEstTimeSeries' object with the estimated signal
    #                covariance
    #     noise.sd:  a vector of the noise standard deviation estimates
    if (!inherits(cov.est.ts, "CovEstTimeSeries")) {
      stop("Argument `cov.est' should be a CovEstTimeSeries object, not: ",
           class(cov.est.ts), ".")
    } else if (!(num.signals.0 >= 0) && num.signals.0 == round(num.signals.0)) {
      stop("Argument `num.signals.0' should be a non-negative integer, not: ",
           num.signals.0, ".")
    } else if (!(noise.sd.0 > 0)) {
      stop("Argument `noise.sd.0' should be positive, not: ",
           noise.sd.0, ".") 
    }

    with(cov.est.ts, {
      rank.est         <- rep(NA, num.times)
      noise.sd.est     <- rep(NA, num.times)
      evalues.sqrt.est <- matrix(NA, num.times, max.rank)
      evectors.est     <- array(NA, c(num.times, dim, max.rank))
     
      for (i in seq_len(num.times)) {
        est <- estimate(num.snapshots[i], dim, num.signals.0, noise.sd.0,
                        evalues[i,], evalues.sqrt[i,],
                        matrix(evectors[i,,], dim, max.rank), ...)
        
        rank.est[i]      <- est$num.signals
        noise.sd.est[i]  <- est$noise.sd
        if (rank.est[i] > 0) {
          evalues.sqrt.est[i,1:rank.est[i]] <- evalues.sqrt[i,1:rank.est[i]]
          evectors.est[i,,1:rank.est[i]]    <- evectors[i,,1:rank.est[i]]
        }

        num.signals.0 <- rank.est[i]
        noise.sd.0    <- noise.sd.est[i]
      }

      cov.est <- MakeCovEstTimeSeries(MakeCovTimeSeries(rank.est,
                                                        evalues.sqrt.est,
                                                        evectors.est,
                                                        time),
                                      num.snapshots)
      list(num.times=num.times, time=time,
           cov.ts=cov.est, noise.sd=noise.sd.est)
    })
  }  
}
