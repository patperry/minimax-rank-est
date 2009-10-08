# landmark-times.R
# ----------------

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

