# cost.R
# ------
# Functions associated with the loss function.
#

source( "rmt.R" )

loss <- function( r, w, lambda, w.est ) {
    b2 <- bias2.est( r, w, lambda, w.est )
    v  <- var.est( r, w, lambda, w.est )
    b2 + v
}

bias2.est <- function( r, w, lambda, w.est ) {
    n <- dim( w )[ 1 ]
    T <- dim( w )[ 3 ]
    
    lambda.sqrt <- sqrt( lambda )
    
    b2 <- matrix( NA, n+1, T )
    for( t in 1:T ) {
        if( r[ t ] > 0 ) {
            w.true <- qr.Q( qr( matrix( w[ ,1:r[ t ],t ], n, r[ t ] ) ) )
            resid  <- w.true %*% Conj(t( w.true ))
            #b2[ 1,t ] <- r[ t ]
            b2[ 1,t ] <- 1
        } else {
            resid <- matrix( 0, n, n )
            b2[ 1,t] <- 0
        }

        for( i in 1:n ) {
            wi <- matrix( w.est[ ,i,t ], n, 1 )
            resid <- resid - wi %*% Conj(t( wi ))
            b2[ i+1,t ] <- sum( sum( abs( resid )^2 ) )
            #b2[ i+1,t ] <- sum( svd( resid )$d )
            #b2[ i+1,t ] <- max( abs( eigen( resid, symmetric=TRUE, only.values=TRUE )$values ) )
        }
    }
    
    b2
}

var.est <- function( r, w, lambda, w.est ) {
    n <- dim( w )[ 1 ]
    T <- dim( w )[ 3 ]
    v <- matrix( 0, n+1, T )
    v
}

# The almost-sure limit of the cost for including an eigenvector in
# the tracked subspace.
#
# Arguments
# ---------
#     lambda : the population eigenvalue
#          c : the aspect ratio ( # of dimensions / sample size )
#     sigma2 : the noise variance
#
# Returns
# -------
#     ci : the cost of including the correpsonding eigenvector
#
cost.include <- function( lambda, c, sigma2=1 )
{
    rbar2 <- eigen.vec.cor2.compl.mean( lambda, c, sigma2 )
    ci    <- ifelse( lambda > 0,
                     rbar2 * lambda + sigma2,
                     sigma2 )
    ci
}

# The almost-sure limit of the cost for excluding an eigenvector from
# the tracked subspace.
#
# Arguments
# ---------
#     lambda : the population eigenvalue
#          c : the aspect ratio ( # of dimensions / sample size )
#     sigma2 : the noise variance
#
# Returns
# -------
#     ce : the cost of excluding the correpsonding eigenvector
#
cost.exclude <- function( lambda, c, sigma2=1 )
{
    ce <- ifelse( lambda > 0, 
                  lambda, 
                  0 )
    ce
}

# The derivative (with respect to the eigenvalue) of the cost for including an
# eigenvector.  See "cost.include".
#
d.cost.include <- function( lambda, c, sigma2=1 )
{
    dci <- ifelse( lambda > sqrt( c )*sigma2,
                   c*(c-1)/(lambda + c)^2,
                   1 )
    dci
}

# The derivitive (with respect to the eigenvalue) of the cost for exluding
# an eigenvector.  See "cost.exclude".
#
d.cost.exclude <- function( lambda, c, sigma2=1 )
{
    dce <- rep( 1, length( lambda ) )
    dce
}

# The value of the population eigenvalue at which the costs for including
# and excluding are equal.
#
# Arguments
# ---------
#          c : the aspect ratio ( # of dimensions / sample size )
#     sigma2 : the noise variance
#
# Returns
# -------
#     lambda : the population eigenvalue at which the costs for including
#             and excluding are equal
#
eigen.val.costs.equal <- function( c, sigma2=1 )
{
    0.5*sigma2*( sqrt( 8*c + 1 ) + 1 )
}

# Min-max optimal rule for deciding whether or not to include a term,
# based on the sample eigenvalue.
#
# Arguments
# ---------
#          l : the sample eigenvalue
#          c : the aspect ratio ( # of dimensions / sample size )
#     sigma2 : the noise variance
#          N : the sample size
#
# Returns
# -------
#     include : TRUE when the rule is to include the term
#
include.eigen.minmax <- function( l, c, sigma2=1, N )
{
    # determine the minmax cutoff point
    lambda.eq <- eigen.val.costs.equal( c, sigma2 )
    dce       <- d.cost.exclude( lambda.eq, c, sigma2 )
    dci       <- d.cost.include( lambda.eq, c, sigma2 )
    ratio     <- min( dce/( dce - dci ), 1 )
        
    eta.eq    <- eigen.val.mean( lambda.eq, c, sigma2 )
    tau2.eq   <- eigen.val.var( lambda.eq, c, sigma2 )
    
    cutoff    <- eta.eq - sqrt( tau2.eq / N ) * qnorm( ratio )

    include   <- l > cutoff
}
