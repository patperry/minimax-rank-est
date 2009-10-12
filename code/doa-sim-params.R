# doa-sim-params.R
# ----------------
source("landmark-times.R")
source("cov-time-series.R")
source("snr-scale.R")

SteeringMatrix <- function(freq.rad, dim) {
  # Get a matrix of steering vectors, assuming a linear array of sensors
  #
  # Args:
  #   freq.rad: the directions of arrival, in radians
  #   dim:      the dimension of the signal (equal to the number of sensors)
  #
  # Returns:
  #   a matrix with the steering vectors as columns (n by r); the column
  #   corresponding to frequency omega has entries
  #     (1/sqrt(n)) * (1, e^(j omega), e^(j 2 omega), ..., e^(j (n-1) omega) ),
  #   where n is the dimension.  Note that this is a unit vector.
  n <- dim
  r <- length(freq.rad)
  k <- 0:(n-1)
  a <- sqrt(1/n) * (matrix(exp(complex(real=0,
                                       imag=kronecker(freq.rad,k))), n, r))
  a
}

MakeDOASimParams <- function(snr.db, freq.deg, dim,
                             time=seq(0, len=num.times),
                             landmarks=kNoLandmarkTimes) {
  # Create a DOA simulation parameter object
  #
  # Args:
  #   snr.db:     a vector of the signal-to-noise ratios (in dB) of the
  #               arrival signals
  #   freq.deg:   a matrix of the arrival angles (in degrees); the i-th row
  #               is a vector of the same length as snr.db giving the arrival
  #               angles of the signals at time i, where `NA' indicates that
  #               a signal is inactive
  #   dim:        the dimension of the simulation
  #   time:       the times of the snapshots (linearly spaced and increasing)
  #   landmarks:  the landmark times (time points when "interesting" activity
  #               occurs) 
  # 
  # Returns:
  #   a `DOASimParams' object initialized with the given values
  num.signals <- length(snr.db)
  num.times   <- nrow(freq.deg)
  if (num.signals != ncol(freq.deg)) {
    stop("Arguments `snr.db' and `freq.deg' have inconsistent dimensions: ",
         length(snr.db), " and ", toString(dim(snr.db)), ".")
  } else if (dim <= 0) {
    stop("Argument `dim' must be positive, not: ", dim, ".")
  } else if (length(time) != num.times) {
    stop("Arguments `freq.deg' and `time' have inconsistent dimensions: ",
         toString(dim(freq.deg)), " and ", length(time), ".") 
  } else if (!inherits(landmarks, "LandmarkTimes")) {
    stop("Argument `landmarks' must be a `LandmarkTimes' object, not: ",
         class(landmarks), ".")
  } 

  freq.rad <- freq.deg *(pi/180)
  snr      <- 10^(snr.db/10)
  snr.sqrt <- 10^(snr.db/20)

  rank         <- rep(NA, num.times)
  evalues.sqrt <- matrix(NA, num.times, num.signals)
  evectors     <- array(NA, c(num.times, dim, num.signals))

  for (i in seq_len(num.times)) {
    active     <- !is.na(freq.rad[i,])
    r          <- sum(active)
    a          <- SteeringMatrix(freq.rad[i,active], dim)
    cov.s.sqrt <- diag(snr.sqrt[active], r)
    
    rank[i] <- r

    if (rank[i] > 0) {
      b                   <- a %*% cov.s.sqrt
      b.svd               <- svd(b, nu=r, nv=0)
      evectors[i,,1:r]    <- b.svd$u
      evalues.sqrt[i,1:r] <- b.svd$d[1:r]
    }
  }
    
  cov.signal <- MakeCovTimeSeries(rank, evalues.sqrt, evectors, time)

  res <- list(num.times=num.times,
                   time=time,
                    dim=dim,
            num.signals=num.signals,
             cov.signal=cov.signal,
                    snr=snr,
               snr.sqrt=snr.sqrt,
                 snr.db=snr.db,
              freq.rad=freq.rad,
              freq.deg=freq.deg,
             landmarks=landmarks)
  
  class(res) <- c("DOASimParams", class(res))
  res 
}


KavcicYang1 <- function(dim=9, freq.snapshot=1) {
  # Get the parameters of the first DOA simulation from Kavcic and Yang (1996),
  # for a given dimension
  #
  # Args:
  #   dim:           the dimension of the simulation
  #   freq.snapthot: the sampling frequency (samples per unit time)
  #
  # Returns:
  #   A `DOASimParams` object with the given simulation parameters.
  #
  # A total of 6 signals appear an disappear in the time interval [0,1000).
  # The signal lifetimes, strengths, and angles of arrival are as follows:
  #
  #    1. [0,300]     0 dB     35 to 25 degrees (linearly interpolated)
  #    2. [0,500]     3 dB     -5 degrees       (constant) 
  #    3. [150,400]  -2 dB     50 degrees       (constant)
  #    4. [150,400]  -2 dB    -30 degrees       (constant)
  #    5. [700,1000]  3 dB     20 to  5 degrees (linearly interpolated)
  #    6. [800,1000]  0 dB    -15 to -5 degrees (linearly interpolated)
  #
  # Samples are given at the integer time points { 0, 1, 2, ..., 999 }.
  if (dim < 6) {
    warning("Argument `dim' is less than the number of signals (6): ",
            dim, ".")
  } else if (!(freq.snapshot > 0)) {
    stop("Argument `freq.snapshot' should be positive, not: ",
         freq.snapshot, ".")
  }

  T     <- 1000
  f     <- freq.snapshot
  index <- function(t) round(t * f) + 1 
  num.snapshots <- index(T)
 
  t0 <-    0; i0 <- index(t0) 
  t1 <-  150; i1 <- index(t1) 
  t2 <-  300; i2 <- index(t2)
  t3 <-  400; i3 <- index(t3)
  t4 <-  500; i4 <- index(t4)
  t5 <-  700; i5 <- index(t5)
  t6 <-  800; i6 <- index(t6)
  t7 <- 1000; i7 <- index(t7)

  snr.db   <- rep(NA, 6)
  freq.deg <- matrix(NA, num.snapshots, 6)
  time     <- seq(0, T, len=num.snapshots)

  freq.deg[i0:i2, 1] <- seq( 35, 25, len=i2-i0+1); snr.db[1] <-  0
  freq.deg[i0:i4, 2] <- rep( -5,     len=i4-i0+1); snr.db[2] <-  3
  freq.deg[i1:i3, 3] <- rep( 50,     len=i3-i1+1); snr.db[3] <- -2
  freq.deg[i1:i3, 4] <- rep(-30,     len=i3-i1+1); snr.db[4] <- -2
  freq.deg[i5:i7, 5] <- seq( 20,  5, len=i7-i5+1); snr.db[5] <-  3
  freq.deg[i6:i7, 6] <- seq(-15, -5, len=i7-i6+1); snr.db[6] <-  0

  landmarks <- MakeLandmarkTimes(c(t0, t1, t2, t3, t4, t5, t6, t7))

  MakeDOASimParams(snr.db, freq.deg, dim, time, landmarks)
}

KavcicYang2 <- function(dim=9, freq.snapshot=1) {
  # Get the parameters of the second DOA simulation from Kavcic and Yang (1996).
  #
  # Args:
  #   dim:           the dimension of the simulation
  #   freq.snapthot: the sampling frequency (samples per unit time)
  #
  # Returns:
  #   A `DOASimParams` object with the given simulation parameters.
  # 
  # A total of 3 signals are present during the time interval [0,1000].
  # The signal lifetimes, strengths, and angles of arrival are as follows:
  # 
  #    1. [0,1000]     4 dB     35 to  20 degrees (linearly interpolated)
  #    2. [0,1000]     4 dB     20 to  35 degrees (linearly interpolated) 
  #    3. [0,1000]    -2 dB      5 to -10 degrees (linearly interpolated)
  #
  if (dim < 3) {
    warning("Argument `dim' is less than the number of signals (3): ",
            dim, ".")
  } else if (!(freq.snapshot > 0)) {
    stop("Argument `freq.snapshot' should be positive, not: ",
         freq.snapshot, ".")
  }
  
  T             <- 1000
  f             <- freq.snapshot
  num.snapshots <- round(T * f) + 1

  time     <- seq(0, T, len=num.snapshots)
  snr.db   <- rep(NA, 3)
  freq.deg <- matrix(NA, num.snapshots, 3)

  freq.deg[,1] <- seq( 35,  20, len=num.snapshots ); snr.db[1] <-  4
  freq.deg[,2] <- seq( 20,  35, len=num.snapshots ); snr.db[2] <-  4
  freq.deg[,3] <- seq(  5, -10, len=num.snapshots ); snr.db[3] <- -2

  landmarks <- MakeLandmarkTimes(c(500))

  MakeDOASimParams(snr.db, freq.deg, dim, time, landmarks)
}

plot.DOASimParams <- function(x, ..., xlab="Time", ylab="Frequency (degrees)",
                              ylim=c(-60,60), lwd=3, 
                              snr.scale=kDefaultSNRScale) {
  with(x, {
    plot(range(time), ylim, t="n", xlab=xlab, ylab=ylab, ...)
    lines(landmarks)
    for (i in seq_len(num.signals)) {
      lines(time, freq.deg[,i], lwd=lwd,
            col=ColorOnSNRScale(snr.scale, snr.db[i]))
    }
  })
}
