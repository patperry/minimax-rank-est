# rank-est-htest.R
# ----------------

require("RMTstat")
source("rank-est-generic.R")

EstimateRankHTest <- function(num.snapshots, dim, num.signals.0, noise.sd.0,
                              evalue, evalue.sqrt, evector, level=0.005,
                              ...) {
  cutoff <- qWishartMax(level, pdim=dim, ndf=num.snapshots, var=noise.sd.0^2,
                        beta=2, lower.tail=FALSE)

  if (num.signals.0 > 0 && evalue[num.signals.0] < cutoff) {
    num.signals <- num.signals.0 - 1
  } else if (num.signals.0 < min(num.snapshots, dim) - 2
             && evalue[num.signals.0+1] > cutoff) {
    num.signals <- num.signals.0 + 1
  } else {
    num.signals <- num.signals.0
  }

  noise.var <- mean(evalue[seq(num.signals + 1, min(dim, num.snapshots))])
  noise.sd  <- sqrt(noise.var)

  list(num.signals=num.signals, noise.sd=noise.sd)
}

EstimateRankHTestTimeSeries <- MakeRankEstTimeSeriesFun(EstimateRankHTest)

