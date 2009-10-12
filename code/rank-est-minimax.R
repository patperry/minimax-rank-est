# rank-est-minimax.R
# ------------------

require("RMTstat")
source("rank-est-generic.R")

MinimaxCutoff <- function(num.snapshots, dim, noise.var=1,
                          spike.signal.min=noise.var*(sqrt(1/svr)+num.snapshots^(-1/3)),
                          cost.exclude=1,
                          cost.include=1,
                          beta=2) {
  svr  <- num.snapshots/dim
  null <- WishartMaxPar(ndf=num.snapshots, pdim=dim, var=noise.var, beta=beta)
  alt  <- WishartSpikePar(spike.signal.min, ndf=num.snapshots, pdim=dim,
                          var=noise.var, beta=beta)

  if (cost.exclude == 0) {
    cutoff <- Inf
  } else if (cost.include == 0) {
    cutoff <- 0
  } else {
    cutoff <- uniroot(interval=null$center + null$scal * c(-10,6), function(l) {
                (cost.include*(1 - ptw((l - null$center)/null$scal, beta=beta))
                 -
                 cost.exclude*pnorm((l - alt$center)/alt$scal))})$root
  }
  cutoff
}

EstimateRankMinimax <- function(num.snapshots, dim, num.signals.0, noise.sd.0,
                                evalue, evalue.sqrt, evector,
                                spike.signal.min=noise.var.0*(sqrt(1/svr)+num.snapshots^(-1/3)),
                                cost.exclude=seq(dim-1, 1, len=dim-1),
                                cost.include=1,
                                beta=2,
                                ...) {
  svr         <- num.snapshots/dim
  noise.var.0 <- noise.sd.0^2

  if (num.signals.0 > 0) {
    cutoff.0 <- MinimaxCutoff(num.snapshots, dim, noise.var.0, spike.signal.min,
                              cost.exclude[num.signals.0], cost.include, beta)
  } else {
    cutoff.0 <- NA
  }

  if (num.signals.0 < min(dim, num.snapshots) - 1) {
    cutoff.1 <- MinimaxCutoff(num.snapshots, dim, noise.var.0, spike.signal.min,
                              cost.exclude[num.signals.0+1], cost.include, beta)
  } else {
    cutoff.1 <- NA
  }

  if (num.signals.0 > 0 && evalue[num.signals.0] < cutoff.0) {
    num.signals <- num.signals.0 - 1
  } else if (num.signals.0 < min(dim, num.snapshots) - 1
             && evalue[num.signals.0+1] > cutoff.1) {
    num.signals <- num.signals.0 + 1
  } else {
    num.signals <- num.signals.0
  }

  noise.var <- mean(evalue[seq(num.signals + 1, min(dim, num.snapshots))])
  noise.sd  <- sqrt(noise.var)

  list(num.signals=num.signals, noise.sd=noise.sd)
}

EstimateRankMinimaxTimeSeries <- MakeRankEstTimeSeriesFun(EstimateRankMinimax)

