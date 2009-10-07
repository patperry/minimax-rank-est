# fit.R
# -----

fit <- function( x, N, sigma2.0, select, ... ) {
    
    n  <- nrow( x )
    T  <- ncol( x )
    r0 <- 0

    cov.window <- function( t ) {
        padding <- rep( 0, max(n-N, 0) )
        
        if( t < N ) {
            x.svd <- svd( x[,1:t], nu=n, nv=0 )
            w     <- x.svd$u
            l     <- c( x.svd$d^2 / t, rep( NA, max( n-t, 0 ) ) )
        } else {
            x.svd <- svd( x[,(t-N+1):t], nu=n, nv=0 )
            w     <- x.svd$u
            l     <- c( x.svd$d^2 / N, padding )
        }
        
        list( w=w, l=l )
    }
    
    r      <- rep( NA, T )
    w      <- array( NA, c(n,n,T) )
    sigma2 <- rep( NA, T )
    l      <- matrix( NA, n, T )
    

    for( t in 1:T ) {
        x.cov <- cov.window( t )

        w[ ,,t ] <- x.cov$w
        l[ ,t  ] <- x.cov$l
        l.signal <- ifelse( r0 > 0,   l[   r0, t ], NA )
        l.noise  <- ifelse( r0 < n-1, l[ r0+1, t ], NA )
        
        r[ t ]   <- select( n=n, N=N, r=r0, sigma2=sigma2.0, 
                        l.signal=l.signal, l.noise=l.noise, ... )
        
        sigma2[ t ] <- ifelse( t < N, sigma2.0, mean( l[ n:(r[ t ]+1), t ] ) )

        r0       <- r[ t ]
        sigma2.0 <- sigma2[ t ]
    }
    
    list( r=r, w=w, sigma2=sigma2, l=l )
}
