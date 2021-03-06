# cov-time-series.R
# -----------------

MakeCovTimeSeries <- function(rank, evalues.sqrt, evectors,
                              time=seq(0, len=num.times)) {
  # Create a covariance matrix time-series object
  #
  # Args:
  #   rank:         a vector of the covariance matrix ranks at each
  #                 time point (num.times)
  #   evalues.sqrt: a matrix of the square roots of the covariance
  #                 eigenvalues at each time point (num.times by max.rank) 
  #   evectors:     an array of the covariance eigenvectors at each
  #                 time point (num.times by dim by max.rank)
  #   time:         a vector of times, in sorted order (num.times)
  #
  # Returns:
  #   a `CovTimeSeries' object initialized with the given values
  num.times <- length(rank)
  max.rank  <- ncol(evalues.sqrt)
  dim       <- dim(evectors)[2] 
  if (num.times != length(time)) {
    stop("Argumets `time' and 'rank' have invalid lenghts: ",
         length(time), " and ", length(rank), ".")
  } else if (num.times != nrow(evalues.sqrt)) {
    stop("Arguments `rank' and `evalues.sqrt' have invalid shapes: ",
         length(rank), " and ", toString(dim(evalues.sqrt)), ".")
  } else if (num.times != dim(evectors)[1]) {
    stop("Arguments `rank' and `evectors' have invalid shapes: ",
         length(rank), " and ", toString(dim(evectors)), ".")
  } else if (max.rank != dim(evectors)[3]) {
    stop("Arguments `evalues.sqrt' and `evectors' have invalid shapes: ",
         toString(dim(evalues.sqrt)), " and ",
         toString(dim(evectors)), ".")
  } else if (!all(0 <= rank & rank <= max.rank & round(rank) == rank)) {
    stop("Argument `rank' should have integer values between 0 and ", max.rank,
         ", not: ", toString(
           rank[which(!(0 <= rank & rank <= max.rank & round(rank) == rank))]),
         ".")
  } else if (!all(evalues.sqrt >= 0 | is.na(evalues.sqrt))) {
    stop("Argument `evalues.sqrt' should have non-negative entries.")
  }

  evalues <- evalues.sqrt^2
  
  res <- list(num.times=num.times,
                    dim=dim,
               max.rank=max.rank,
                   time=time,
                   rank=rank,
                evalues=evalues,
           evalues.sqrt=evalues.sqrt,
               evectors=evectors)
  class(res) <- c("CovTimeSeries", class(res))
  res              
}

MakeCovEstTimeSeries <- function(cov.ts, num.snapshots) {
  # Annotate an estimated "CovTimeSeries" object with the number of snapshots
  # used to compute the estimates.
  # 
  # Args:
  #   cov.ts:        a CovTimeSeries ojbect
  #   num.snapshots: a vector of the number of snapshots
  #
  # Returns:
  #   a `CovEstTimeSeries' object.
  if (!inherits(cov.ts, "CovTimeSeries")) {
    stop("Argument `cov.ts' should be a CovTimeSeries object, not: ",
         class(cov.ts), ".")
  } else if (cov.ts$num.times != length(num.snapshots)) {
    stop("Arguments `cov.ts' and `num.snapshots' have inconsistent lengths: ",
         cov.ts$num.times, " and ", length(num.snapshots), ".")
  }

  res <- c(cov.ts, list(num.snapshots=num.snapshots))
  class(res) <- c("CovEstTimeSeries", class(cov.ts))
  res
}

SampleCovTimeSeries <- function(cov.ts) {
  # Sample gaussian snapshots with the given covariances
  #
  # Args:
  #   cov.ts: a CovTimeSeries object
  #
  # Returns:
  #   a matrix whose rows are Gaussian obervations; the covariance of the
  #   `i`th row is gotten from the `i`th timepoint in `cov.ts`.
  if (!inherits(cov.ts, "CovTimeSeries")) {
    stop("Argument `cov.ts' should be a CovTimeSeries object, not: ",
         class(cov.ts), ".")
  }

  num.snapshots <- cov.ts$num.times
  dim           <- cov.ts$dim
  x             <- matrix(NA, num.snapshots, dim)

  if (is.complex(cov.ts$evectors)) {
    rgauss <- function(n, sd) complex(real=rnorm(n, sd=sd*sqrt(1/2)),
                                      imag=rnorm(n, sd=sd*sqrt(1/2)))
  } else {
    rgauss <- rnorm
  }

  for (i in seq_len(num.snapshots)) {
    r <- cov.ts$rank[i]
    if (r > 0) {
      d     <- cov.ts$evalues.sqrt[i,1:r]
      w     <- matrix(cov.ts$evectors[i,,1:r], dim, r) 
      z     <- rgauss(r, sd=d)
      x[i,] <- w %*% cbind(z) 
    } else {
      x[i,] <- 0
    }
  }
  
  x
}


WindowedCovEst <- function(snapshot, length.window, 
                           time=seq(0, len=num.snapshots)) {
  # Compute a windowed covariance estimate from centered (mean-zero) data
  #
  # Args:
  #   snapshot:      a matrix of snapshots (num.snapshots by dim)
  #   length.window: the length of the window, in snapshots
  #   time:          times associated with the snapshots (num.snapshots)
  #
  # Returns:
  #   a `CovTimeSeries' object with the windowed estimate 
  snapshot      <- as.matrix(snapshot)
  num.snapshots <- nrow(snapshot)
  dim           <- ncol(snapshot)
  max.rank      <- min(dim, length.window)

  num.window.snapshots <- rep(NA, num.snapshots)
  rank                 <- rep(NA, num.snapshots)
  evalues.sqrt         <- matrix(NA, num.snapshots, max.rank)
  evectors             <- array(NA, c(num.snapshots, dim, max.rank))
 
  for (i in seq_len(min(length.window, num.snapshots))) {
    num.window.snapshots[i] <- i
    r                       <- min(i,max.rank)
    x                       <- snapshot[1:i,,drop=FALSE]
    x.svd                   <- svd(x, nu=0, nv=r)

    rank[i]                 <- r
    evalues.sqrt[i,1:r]     <- x.svd$d[1:r] / sqrt(i)
    evectors[i,,1:r]        <- x.svd$v
  }

  r <- min(length.window, dim)
  for (i in seq(length.window+1, len=max(num.snapshots - length.window, 0))) {
    num.window.snapshots[i] <- length.window
    x                       <- snapshot[(i-length.window+1):i,]
    x.svd                   <- svd(x, nu=0, nv=r)  
    
    rank[i]          <- r
    evalues.sqrt[i,] <- x.svd$d[1:r] / sqrt(length.window)
    evectors[i,,]    <- x.svd$v
  }

  cov.ts <- MakeCovTimeSeries(rank, evalues.sqrt, evectors, time)
  MakeCovEstTimeSeries(cov.ts, num.window.snapshots)
}

DiffCovTimeSeriesSubspaceFrob2 <- function(cov.ts1, cov.ts2) {
  if (!inherits(cov.ts1, "CovTimeSeries")) {
    stop("Argument `cov.ts1' should be a `CovTimeSeries' object, not: ",
         toString(class(cov.ts1)), ".")
  } else if (!inherits(cov.ts2, "CovTimeSeries")) {
    stop("Argument `cov.ts2' should be a `CovTimeSeries' object, not: ",
         toString(class(cov.ts2)), ".")
  } else if (!identical(cov.ts1$time, cov.ts2$time)) {
    stop("Arguments `cov.ts1' and `cov.ts2' should have the same time points.")
  } else if (cov.ts1$dim != cov.ts2$dim) {
    stop("Arguments `cov.ts' and `cov.ts2' have inconsistent dimensions: ",
         cov.ts1$dim, " and ", cov.ts2$dim, ".")
  }

  num.times <- cov.ts1$num.times
  time      <- cov.ts1$time
  dim       <- cov.ts1$dim
  frob2     <- rep(NA, num.times)

  for (i in seq_len(num.times)) {
    r1 <- cov.ts1$rank[i]
    if (r1 > 0) {
      w1 <- cov.ts1$evectors[i,,1:r1]
      p1 <- w1 %*% Conj(t(w1))
    } else {
      p1 <- matrix(0, dim, dim)
    }
    
    r2 <- cov.ts2$rank[i]
    if (r2 > 0) {
      w2 <- cov.ts2$evectors[i,,1:r2]
      p2 <- w2 %*% Conj(t(w2))
    } else {
      p2 <- matrix(0, dim, dim)
    }

    delta <- p1 - p2
    frob2[i] <- sum(abs(delta)^2)
  }

  list(num.times=num.times, time=time, frob2=frob2)
}
