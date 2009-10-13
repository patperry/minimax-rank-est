# rank-est-oracle.R
# -----------------

source("cov-time-series.R")

EstimateRankOracleTimeSeries <- function(cov.true, cov.est) {
  if (!inherits(cov.true, "CovTimeSeries")) {
    stop("Argument `cov.true' should be a `CovTimeSeries object, not: ",
         class(cov.true), ".")
  } else if (!inherits(cov.est, "CovEstTimeSeries")) {
    stop("Argument `cov.est' should be a `CovEstTimeSeries object, not: ",
         class(cov.est), ".")
  } else if (!identical(cov.true$time, cov.est$time)) {
    stop("Arguments `cov.true' and `cov.est' should have the same time pionts.")
  } else if (cov.true$dim != cov.est$dim) {
    stop("Arguments `cov.true' and `cov.est' should have the same `dim', not: ",
         cov.true$dim, " and ", cov.est$dim, ".")
  }

  num.times     <- cov.true$num.times
  dim           <- cov.true$dim
  max.rank      <- cov.true$max.rank
  time          <- cov.true$time
  num.snapshots <- cov.est$num.snapshots

  rank         <- rep(NA, num.times)
  evalues.sqrt <- matrix(NA, num.times, max.rank)
  evectors     <- array(NA, c(num.times, dim, max.rank))
  noise.sd     <- rep(NA, num.times)

  for (i in seq_len(num.times)) {
    r       <- cov.true$rank[i]
    rank[i] <- r

    r <- min(r, num.snapshots)
    if (r > 0) {
      evalues.sqrt[i,1:r] <- cov.est$evalues.sqrt[i,1:r]
      evectors[i,,1:r]    <- cov.est$evectors[i,,1:r]
    } 
    if (r < rank[i]) {
      evalues.sqrt[i,(r+1):rank[i]] <- 0
      evectors[i,,(r+1):rank[i]]    <- 0
    }
    
    if (r < dim) {
      noise.sd[i] <- sqrt(mean(cov.est$evalues[i,(r+1):dim]))
    } else {
      noise.sd[i] <- 0
    }
  }

  cov.ts <- MakeCovEstTimeSeries(MakeCovTimeSeries(rank,
                                                   evalues.sqrt,
                                                   evectors,
                                                   time),
                                 num.snapshots)
  list(num.times=num.times, time=time, cov.ts=cov.ts, noise.sd=noise.sd) 
}
