# doa-sim-params.R
# ----------------
source("landmark-times.R")

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
  # Get the active signal strengths and directions at particular times
  #
  # Args:
  #   par.sim: a DOASimParams object with simulation parameters
  #   t:       a vector of times
  #
  # Returns:
  #   a list with the following components:
  #     active:   a boolean vector indicating which signals are active
  #     snr:      the strenghs (signal-to-noise ratios) of all signals
  #     snr.sqrt: the square-roots of the strengths
  #     snr.db:   the strengths in decibels
  #     freq.deg: the directions of arrival of the signals (in degrees)
  #     freq.rad: the directions of arrival in radians
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

