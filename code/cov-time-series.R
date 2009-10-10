# cov-time-series.R
# -----------------

MakeCovTimeSeries <- function(time, rank, evalues.sqrt, evectors) {
  # Create a covariance matrix time-series object
  #
  # Args:
  #   time:         a vector of times, in sorted order (N)
  #   rank:         a vector of the covariance matrix ranks at each
  #                 time point (N)
  #   evalues.sqrt: a matrix of the square roots of the covariance
  #                 eigenvalues at each time point (N by max.rank) 
  #   evectors:     an array of the covariance eigenvectors at each
  #                 time point (N by n by max.rank)
  #
  # Returns:
  #   a `CovTimeSeries' object initialized with the given values
  num.times <- length(time)
  max.rank  <- ncol(evalues.sqrt)
  dim       <- dim(evectors)[2] 
  if (num.times != length(rank)) {
    stop("Argumets `time' and 'rank' have invalid lenghts: ",
         length(time), " and ", length(rank), ".")
  } else if (num.times != nrow(evalues.sqrt)) {
    stop("Arguments `time' and `evalues.sqrt' have invalid shapes: ",
         length(time), " and ", toString(dim(evalues.sqrt)), ".")
  } else if (num.times != dim(evectors)[1]) {
    stop("Arguments `time' and `evectors' have invalid shapes: ",
         length(time), " and ", toString(dim(evectors)), ".")
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

WindowedCovEst <- function(snapshot, length.window, 
                           time=seq(0, len=num.snapshots)) {
  # Compute a windowed covariance estimate from centered (mean-zero) data
  #
  # Args:
  #   snapshot:      a matrix of snapshots (N by n)
  #   length.window: the length of the window
  #   time:          times associated with the snapshots (N)
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

  MakeCovTimeSeries(time, rank, evalues.sqrt, evectors)
}
