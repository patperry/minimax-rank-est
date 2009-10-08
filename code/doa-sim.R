# doa-sim.R
# ---------
source("doa-sim-params.R")
source("cov-time-series.R")


SampleDOASim <- (function() {
  steering <- function(freq.rad, dim) {
    n <- dim
    r <- length(freq.rad)
    k <- 0:(n-1)
    a <- sqrt(1/n) * (matrix(exp(complex(real=0,
                                         imag=kronecker(freq.rad,k))), n, r))
    a
  }

  function(par.sim, dim, freq.snapshot=1) {
    # Sample a DOA simulation object
    #
    # Args:
    #   par.sim:       a DOASimParams object with the simulation parameters
    #   dim:           the dimension of the simulation
    #   freq.snapshot: the snapshot sampling frequency (snapshots per unit
    #                  time)
    #
    # Returns:
    #   a `DOASim' object initialized with the given values
    if (!inherits(par.sim, "DOASimParams")) {
      stop("Argument `par.sim' must be a `DOASimParams' object, not: ",
           class(par.sim), ".") 
    } else if (!(round(dim) == dim && dim > 0)) { 
      stop("Argument `dim' must be a positive integer, not: ", dim, ".")
    } else if (!(freq.snapshot > 0)) {
      stop("Argument `freq.snapshot' must be positive, not: ",
           freq.snapshot, ".")
    }

    T           <- par.sim$max.time
    num.signals <- par.sim$num.signals
    n           <- dim
    f           <- freq.snapshot
    Tf          <- round(T * f) + 1

    time   <- seq(0, T, len=Tf)
    x      <- matrix(NA, Tf, n)
    rank   <- rep(NA, Tf)
    w      <- array(NA, c(Tf, n, num.signals))
    lambda.sqrt <- matrix(NA, Tf, num.signals)

    e <- matrix(complex(real=rnorm(Tf*n, sd=sqrt(1/2)),
                        imag=rnorm(Tf*n, sd=sqrt(1/2))),
                Tf, n)
    s <- matrix(complex(real=rnorm(Tf*num.signals, sd=sqrt(1/2)),
                        imag=rnorm(Tf*num.signals, sd=sqrt(1/2))),
                Tf, num.signals)

    for (i in seq_len(Tf)) {
      t       <- time[i]
      doa     <- DOAAtTime(par.sim, t)
      r       <- sum(doa$active)
      a       <- (steering(doa$freq.rad[doa$active], n) 
                  %*% diag(doa$snr.sqrt[doa$active], r))
      rank[i] <- r

      if (rank[i] > 0) {
        a.svd              <- svd(a, nu=r, nv=r)
        w[i,,seq_len(r)]   <- a.svd$u
        lambda.sqrt[i,1:r] <- a.svd$d[1:r]
        x[i,]              <- (a.svd$u
                               %*% (diag(lambda.sqrt[i,1:r], r)
                                    %*% cbind(s[i,doa$active][1:r]))
                               + e[i,])
      } else {
        x[i,] <- e[i,]
      }
 
      s[i,!doa$active] <- NA
    }

    cov.signal <- MakeCovTimeSeries(time,rank,lambda.sqrt,w)

    res <- list(par=par.sim,
                dim=n,
      freq.snapshot=freq.snapshot,
      num.snapshots=Tf,
               time=time,
         cov.signal=cov.signal,
           snapshot=x,
     signal.sphered=s,
              noise=e)
    
    class(res) <- c("DOASim", class(res))
    res
  }
})()

