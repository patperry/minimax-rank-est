# select.R
# --------
# Decision rules for selecting the rank of the tracked subspace.
#

source( "cost.R" )
require( "RMTstat", "~/Lib/R" )

# Esitmate the rank of the subspace, using a min-max optimal decision rule.  
# Note that the rule is only valid when $n < m$.
# 
# Arguments
# ---------
#            n : the dimensionality of the data
#            N : the window size
#           r0 : the current rank
#       sigma2 : the noise variance
#     l.signal : the smallest signal eigenvalue
#      l.noise : the largest noise eigenvalue
# 
# Returns
# -------
#     r : the rank of the subspace, one of { r0-1, r0, r0+1 }.
#
select.rank.minmax <- function( n, N, r0, sigma2, l.signal, l.noise ) {
    c <- n/N
    
    if( r0 > 0 && !include.eigen.minmax( l.signal, c, sigma2, N ) ) {
        r <- r0 - 1
    } else if( r0 < n-1 && include.eigen.minmax( l.noise, c, sigma2, N ) ) {
        r <- r0 + 1
    } else {
        r <- r0
    }
    
    r
}

# Esitmate the rank of the subspace, using plugin estimates for the costs
# of including or excluding terms.
# 
# Arguments
# ---------
#            n : the dimensionality of the data
#            N : the window size
#           r0 : the current rank
#       sigma2 : the noise variance
#     l.signal : the smallest signal eigenvalue
#      l.noise : the largest noise eigenvalue
# 
# Returns
# -------
#     r : the rank of the subspace, one of { r0-1, r0, r0+1 }.
#
select.rank.plain <- function( n, N, r0, sigma2, l.signal, l.noise ) {
    c      <- n/N
    k      <- eigen.val.costs.equal( c, sigma2 )
    cutoff <- eigen.val.mean( k, c, sigma2 )
    
    if( r0 > 0 && l.signal < cutoff ) {
        r <- r0 - 1
    } else if( r0 < n-1 && l.noise > cutoff ) {
        r <- r0 + 1
    } else {
        r <- r0
    }
    
    r
}

select.tw <- function( n, N, r0, sigma2, l.signal, l.noise ) {
    gamma  <- n/N
    params <- wishart.max.par( N, n, sigma2, beta=2 )
    t      <- qtw( 0.5, beta=2 )
    cutoff <- params$center + t*params$scale

    if( r0 > 0 && l.signal < cutoff ) {
        r <- r0 - 1
    } else if( r0 < n-1 && l.noise > cutoff ) {
        r <- r0 + 1
    } else {
        r <- r0
    }
    
    r    
}

minmax.cutoff <- function( n, N, sigma2, lambda.min=10^(-2/10), 
                           cost.include=1, cost.exclude=1 ) {
    gamma        <- n/N
    lambda.min   <- 10^(-2/10)
    
    tw    <- wishart.max.par( N, n, sigma2, beta=2 )
    spike <- wishart.spike.par( lambda.min, N, n, sigma2 )
    spike$scale <- spike$scale/sqrt(2)
    
    T <- uniroot( interval=c( -10, 6 ),
             function( t ) { 
                 ( ( cost.include
                     * (1 - ptw( ( t - tw$center )/tw$scale, beta=2 ) ) )
                   -
                   ( cost.exclude
                     * pnorm( ( t - spike$center )/spike$scale ) ) )
             })$root
    T
}

minmax.cutoffs <- function( n, N, sigma2, lambda.min=(sqrt(n/N)+N^(-1/3)), 
                           cost.include=1, cost.exclude=1 ) {
    T <- rep( NA, n )
    for( i in 1:n ) {
        T[ i ] <- minmax.cutoff( n, N, sigma2, lambda.min, 
                                 cost.include, cost.exclude*(n-i+1) )
    }
    
    T
}


select.minmax <- function( n, N, r0, sigma2, l.signal, l.noise ) {
    cutoffs <- minmax.cutoffs( n, N, sigma2 )

    if( r0 > 0 && l.signal < cutoffs[ r0 ] ) {
        r <- r0 - 1
    } else if( r0 < n-1 && l.noise > cutoffs[ r0+1 ] ) {
        r <- r0 + 1
    } else {
        r <- r0
    }
    
    r    
}

select.kn <- function( n, N, r0, sigma2, l.signal, l.noise ) {
    T <- qwishart.max(0.005, p.dim=n, n.df=N, var=sigma2, beta=2, lower.tail=FALSE)
    
    if( r0 > 0 && l.signal < T ) {
        r <- r0 - 1
    } else if( r0 < n-1 && l.noise > T ) {
        r <- r0 + 1
    } else {
        r <- r0
    }
    
    r    
}

select.50 <- function( n, N, r0, sigma2, l.signal, l.noise ) {
    T <- qwishart.max(0.5, p.dim=n, n.df=N, var=sigma2, beta=2)
    
    if( r0 > 0 && l.signal < T ) {
        r <- r0 - 1
    } else if( r0 < n-1 && l.noise > T ) {
        r <- r0 + 1
    } else {
        r <- r0
    }
    
    r    
}

# Default selection rule.
#select.rank <- select.tw
#select.rank <- select.minmax
select.rank <- select.kn

