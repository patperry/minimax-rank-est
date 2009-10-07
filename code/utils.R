# utils.R
# -------
# Utility functions.

# Generate a matrix of steering vectors.
#
# Arguments
# ---------
#          dim     : the dimensionality of the data (n)
#          snr (r) : a vector of signal-to-noise ratios (*not* in decibels)
#     freq.rad (r) : a vector a steering frequencies, in radians.  In frequency
#                   retrieval, @freq.rad[i]@ is the angular frequency of
#                   the @i@th sinusoid.  In array processing, 
#                       @freq.rad[i] = 2 pi (d/lambda) sin( theta[i] )@,
#                   holds when plane waves impinge a linear uniform sensor
#                   array.  Here @d@ is the spacing between adjacent sensor
#                   elements, @lambda@ is the wavelength, and @theta[i]@ is
#                   the DOA relative to the array broadside.
#
# Returns
# -------
#     a (n*r) : a matrix whose @i@th column contains the @i@th steering vector,
#              scaled by the signal-to-noise ratio.  The steering vector for
#              frequency @omega@ is given by 
#                  @(1, exp(j*omega), ..., exp(j*(n-1)*omega))/sqrt(n)@.
#
steering <- function( dim, snr, freq.rad ) {
    n     <- dim
    nfreq <- length( freq.rad )
    
    k <- 0:(n-1)
    a <- matrix( exp( complex( real=0, imag=kronecker( freq.rad, k ) ) ),
                 n, nfreq )
    a %*% diag( sqrt( snr/n ), nfreq, nfreq )
}
