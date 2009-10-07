# init.R
# ------
# Compute initial values for a subspace tracker.
#

source( "utils.R" )

# Generate appropriate initial values for tracking a time-varying
# low-rank + noise model.  
#
# Arguments
# ---------
#     dim                 : the dimension of the snapshots (n)
#     freq.deg (r.total)  : a vector of steering frequencies, the $t$th
#                           column giving the frequencies of the different
#                           components (in degrees).  'NA' indicates
#                           inactive components.
#     snr      (r.total)  : a vector of component strengths, in decibels
#
# Returns
# -------
#     r0                  : the initial number of components
#     w0 (n * r0)         : an orthonormal basis for the components
#     sigma2.0            : the initial noise level
#     c0 (n * n)          : the initial snapshot covariance
#
init.subspace <- function( dim, freq.deg, snr.db )
{
    n <- dim
    
    missing  <- is.na( freq.deg ) | is.na( snr.db )
    freq.deg <- freq.deg[ !missing ]
    snr.db   <- snr.db[ !missing ]
    
    freq.rad <- freq.deg*( pi/180 )
    snr      <- 10^( snr.db/10 )
    
    r0       <- sum( !missing )
    w0       <- qr.Q( qr( steering( n, snr=1, freq.rad ) ) )
    sigma2.0 <- 1
    c0       <- ( w0 %*% diag( snr, r0, r0 ) %*% Conj( t( w0 ) ) 
                  + diag( sigma2.0, n, n ) )
    
    list( r0=r0, w0=w0, sigma2.0=sigma2.0, c0=c0 )
}
