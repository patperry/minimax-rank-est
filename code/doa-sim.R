
MakeLandmarkTimes <- function(t.landmark) {
  # Create a LandmarkTimes object
  #
  # Args:
  #   t.landmark: a vector giving the times of the landmarks
  #
  # Returns:
  #   a `LandmarkTimes' object initialized with the given values

  if (length(t.landmark) > 0) {
    res <- unique(sort(t.landmark))
  } else {
    res <- c()
  }
  res <- as.numeric(res)

  class(res) <- c("LandmarkTimes", class(res))
  res 
}

# Constant: kNoLandmarkTimes
#   A LandmarkTimes object with no landmarks.
kNoLandmarkTimes <- MakeLandmarkTimes(c())


MakeDOASimParams <- function(snr.db, freq.deg, t.landmark=kNoLandmarkTimes) {
  # Create a DOA simulation parameter object
  #
  # Args:
  #   snr.db:     a vector of the signal-to-noise ratios (in dB) of the
  #               arrival signals
  #   freq.deg:   a matrix of the arrival angles (in degrees); the i-th row
  #               is a vector of the same length as snr.db giving the arrival
  #               angles of the signals at time i, where `NA' indicates that
  #               a signal is inactive
  #   t.landmark: the landmark times (time points when "interesting" activity
  #               occurs) 
  # 
  # Returns:
  #   a `DOASimParams' object initialized with the given values
  r.max <- length(snr.db)
  if (r.max != ncol(freq.deg)) {
    stop("Arguments `snr.db' and `freq.deg' have inconsistent dimensions: ",
         length(snr.db), "and ", dim(snr.db), ".")
  }
  if (!inherits(t.landmark, "LandmarkTimes")) {
    stop("Argument `t.landmark' must be a `LandmarkTimes' object, not: ",
         class(t.landmark), ".")
  } 
  
  T <- nrow(freq.deg) - 1

  freq.rad <- freq.deg *(pi/180)
  snr      <- 10^(snr.db/10)
  snr.sqrt <- 10^(snr.db/20)

  res <- list(max.time=T,
           num.signals=r.max,
                   snr=snr,
              snr.sqrt=snr.sqrt,
                snr.db=snr.db,
              freq.rad=freq.rad,
              freq.deg=freq.deg,
            t.landmark=t.landmark)
  class(res) <- c("DOASimParams", class(res))
  res 
}

DOAAtTime <- function(par.sim, t) {
  if (!inherits(par.sim, "DOASimParams")) {
    stop("Argument `par.sim' should be a DOASimParams object, not: ",
         class(par.sim), ".")
  } 
  
  T <- par.sim$max.time
  if (!all(0 <= t & t <= T)) {
    stop("Error in argument `t'; ",
         "times must be in the range [0,", T, "], not: ",
         toString(t[which(!(0 <= t & t <= T))]), ".")
  }

  num.times   <- length(t)
  num.signals <- par.sim$num.signals
  snr      <- drop(matrix(par.sim$snr,      byrow=TRUE, num.times, num.signals))
  snr.sqrt <- drop(matrix(par.sim$snr.sqrt, byrow=TRUE, num.times, num.signals))
  snr.db   <- drop(matrix(par.sim$snr.db,   byrow=TRUE, num.times, num.signals))
  
  t.prev <- floor(t)
  t.next <- floor(t+1)
  w.prev <- (t.next - t)
  w.next <- (t - t.prev)
  t.next <- pmin(t.next, T) # make sure `interp' works for t=T
  interp <- function(x) w.prev * x[t.prev+1,] + w.next * x[t.next+1,]

  freq.rad <- interp(par.sim$freq.rad)
  freq.deg <- interp(par.sim$freq.deg)

  active <- !is.na(freq.rad)

  res <- list(active=active,
                 snr=snr,
            snr.sqrt=snr.sqrt,
              snr.db=snr.db,
            freq.rad=freq.rad,
            freq.deg=freq.deg)
  res 
}


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
    lambda <- matrix(NA, Tf, num.signals)

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
        a.svd                <- svd(a, nu=r, nv=r)
        a.svd$d              <- c(a.svd$d, rep(0, max(0, r-n)))
        w[i,,seq_len(r)]     <- a.svd$u
        lambda.sqrt          <- a.svd$d[1:rank[i]]
        lambda[i,seq_len(r)] <- lambda.sqrt^2
        x[i,]                <- (a.svd$u
                                 %*% (diag(lambda.sqrt, r)
                                      %*% cbind(s[i,doa$active]))
                                 + e[i,])
      } else {
        x[i,] <- e[i,]
      }
 
      s[i,!doa$active]     <- NA
    }

    cov.signal <- list(evectors=w, evalues=lambda)

    res <- list(par=par.sim,
                dim=n,
      freq.snapshot=freq.snapshot,
      num.snapshots=Tf,
               time=time,
           snapshot=x,
         cov.signal=cov.signal,
               rank=rank,
     signal.sphered=s,
              noise=e)
    
    class(res) <- c("DOASim", class(res))
    res
  }
})()


Sample <- function(object, ...) UseMethod("Sample")
Sample.DOASimParams <- function(object, dim, freq=1, ...) {
  # Generate snapshots 
}


# Constant: kKavcicYang1
#   The parameters of the first DOA simulation from Kavcic and Yang (1996).
#
#   A total of 6 signals appear an disappear in the time interval [0,1000).
#   The signal lifetimes, strengths, and angles of arrival are as follows:
#
#      1. [0,300]     0 dB     35 to 25 degrees (linearly interpolated)
#      2. [0,500]     3 dB     -5 degrees       (constant) 
#      3. [150,400]  -2 dB     50 degrees       (constant)
#      4. [150,400]  -2 dB    -30 degrees       (constant)
#      5. [700,1000]  3 dB     20 to  5 degrees (linearly interpolated)
#      6. [800,1000]  0 dB    -15 to -5 degrees (linearly interpolated)
#
#   Samples are given at the integer time points { 0, 1, 2, ..., 999 }.
kKavcicYang1 <- (function() {
  T  <- 1000

  t0 <-    0
  t1 <-  150    
  t2 <-  300
  t3 <-  400
  t4 <-  500
  t5 <-  700
  t6 <-  800
  t7 <- 1000

  snr.db   <- rep(NA, 6)
  freq.deg <- matrix(NA, T+1, 6)

  freq.deg[1+t0:t2, 1] <- seq( 35, 25, len=t2-t0+1); snr.db[1] <-  0
  freq.deg[1+t0:t4, 2] <- rep( -5,     len=t4-t0+1); snr.db[2] <-  3
  freq.deg[1+t1:t3, 3] <- rep( 50,     len=t3-t1+1); snr.db[3] <- -2
  freq.deg[1+t1:t3, 4] <- rep(-30,     len=t3-t1+1); snr.db[4] <- -2
  freq.deg[1+t5:t7, 5] <- seq( 20,  5, len=t7-t5+1); snr.db[5] <-  3
  freq.deg[1+t6:t7, 6] <- seq(-15, -5, len=t7-t6+1); snr.db[6] <-  0

  t.landmark <- MakeLandmarkTimes(c(t0, t1, t2, t3, t4, t5, t6, t7))

  MakeDOASimParams(snr.db, freq.deg, t.landmark)
})()

# Constant: kKavcicYang2
#   The parameters of the second DOA simulation from Kavcic and Yang (1996).
#
#   A total of 3 signals are present during the time interval [0,1000).
#   The signal lifetimes, strengths, and angles of arrival are as follows:
#
#      1. [0,1000]     4 dB     35 to  20 degrees (linearly interpolated)
#      2. [0,1000]     4 dB     20 to  35 degrees (linearly interpolated) 
#      3. [0,1000]    -2 dB      5 to -10 degrees (linearly interpolated)
#
#   Samples are given at the integer time points { 0, 1, 2, ..., 999 }.
kKavcicYang2 <- (function() {
  T <- 1000

  snr.db   <- rep(NA, 3)
  freq.deg <- matrix(NA, T+1, 3)

  freq.deg[,1] <- seq( 35,  20, len=T+1 ); snr.db[1] <-  4
  freq.deg[,2] <- seq( 20,  35, len=T+1 ); snr.db[2] <-  4
  freq.deg[,3] <- seq(  5, -10, len=T+1 ); snr.db[3] <- -2

  t.landmark <- MakeLandmarkTimes(c(500))

  MakeDOASimParams(snr.db, freq.deg, t.landmark)
})()

