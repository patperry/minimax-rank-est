# generate.R
# ----------
# Generate a data stream from a time-varying low-rank + noise model.
#

source( "utils.R" )

# Generate a sequence of $n$-dimensional snapshots from a time-varying
# low-rank + noise model.  There are at most $r.total$ different components
# active active over the course of the model.  When active, each component
# has a constant  amplitude, but its steering frequency can be time-varying.
#
# Arguments
# ---------
#     dim                    : the dimension of the snapshots (n)
#     freq.deg (r.total * T) : a matrix of steering frequencies, the $t$th
#                             column giving the frequencies of the different
#                             components (in degrees).  'NA' indicates
#                             inactive components.
#    
#     snr      (r.total)     : a vector of component strengths, in decibels
#
# Returns
# -------
#     x (n * T)       : the matrix of snapshots
#     s (n * r.total) : the matrix of steering weights
#     e (n * T)       : the matrix of background noise
#
generate <- function( dim, freq.deg, snr.db ) {
    T       <- ncol( freq.deg )
    r.total <- nrow( freq.deg )
    n       <- dim
    
    freq.rad                      <- freq.deg*( pi/180 )
    snr                           <- 10^( snr.db/10 )
    
    e <- matrix( complex( real=rnorm( n*T, sd=sqrt( 1/2 ) ),
                          imag=rnorm( n*T, sd=sqrt( 1/2 ) ) ), n, T )
    s <- matrix( complex( real=rnorm( r.total*T, sd=sqrt( 1/2 ) ),
                          imag=rnorm( r.total*T, sd=sqrt( 1/2 ) ) ), r.total, T )

    x      <- matrix( NA, n, T )
    r      <- rep( NA, T )
    w      <- array( NA, c(n, r.total, T) )
    lambda <- matrix( NA, r.total, T )
    
    for( t in 1:T ) {
        active  <- !is.na( freq.rad[,t] )
        
        r[ t ]  <- sum( active )
        if( r[ t ] > 0 ) {
            a                  <- steering( n, 1, freq.rad[ active,t ] )
            ul                 <- a %*% diag( sqrt( snr[ active ] ), r[ t ] )
            ul.svd             <- svd( ul )
            w[ ,1:r[t],t ]     <- ul.svd$u[ ,1:r[t] ]
            lambda[ 1:r[t],t ] <- ( ul.svd$d[ 1:r[t] ] )^2
            x[ ,t ]            <- ul %*% s[ active,t ] + e[ ,t ]
        } else {
            x[ ,t ] <- e[ ,t ]
        }
    }
    
    list( x=x, r=r, w=w, lambda=lambda, s=s, e=e )
}
