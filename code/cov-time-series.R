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

  rank         <- rep(NA, num.snapshots)
  evalues.sqrt <- matrix(NA, num.snapshots, max.rank)
  evectors     <- array(NA, c(num.snapshots, dim, max.rank))
 
  for (i in seq_len(min(length.window, num.snapshots))) {
    r       <- min(i,max.rank)
    x       <- snapshot[1:i,,drop=FALSE]
    x.svd   <- svd(x, nu=0, nv=r)

    rank[i]             <- r
    evalues.sqrt[i,1:r] <- x.svd$d[1:r] / sqrt(i)
    evectors[i,,1:r]    <- x.svd$v
  }

  r <- min(length.window, dim)
  for (i in seq(length.window+1, len=max(num.snapshots - length.window, 0))) {
    x     <- snapshot[(i-length.window+1):i,]
    x.svd <- svd(x, nu=0, nv=r)  
    
    rank[i]          <- r
    evalues.sqrt[i,] <- x.svd$d[1:r] / sqrt(length.window)
    evectors[i,,]    <- x.svd$v
  }

  MakeCovTimeSeries(rank, evalues.sqrt, evectors, time)
}

