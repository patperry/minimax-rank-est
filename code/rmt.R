# rmt.R
# -----
# Quantities related to the eigenvalues and eigenvectors from sample
# covariance matrices in spiked models.  The true covariance is assumed
# to be of the form $\Sigma = U \Lambda U^* + \sigma^2 I$, where $U$ is
# a low-rank matrix.
#


# Almost-sure limit of a sample eigenvalue.
#
# Arguments
# ---------
#     lambda : the population eigenvalue
#          c : the aspect ratio ( # of variables / sample size )
#     sigma2 : the noise variance
#
# Returns
# -------
#     eta : the almost-sure limit of the sample eigenvalue
#
eigen.val.mean <- function( lambda, c, sigma2=1 )
{
    eta <- ifelse( lambda > sqrt( c )*sigma2,
                   (lambda + sigma2)*(1 + c*sigma2/lambda),
                   (sqrt( c ) + 1)^2 * sigma2 )
    eta
}

# Scaled limit of the variance of a sample eigenvalue.  
#
# Arguments
# ---------
#     lambda : the population eigenvalue
#          c : the aspect ratio ( # of variables / sample size )
#     sigma2 : the noise variance
#
# Returns
# -------
#     tau2 : the scaled variance of the sample eigenvalue.  This value
#           divided by sample size is equal to the variance of the sample
#           eigenvalue (with error o(1/n), where n is the sample size).
#
eigen.val.var <- function( lambda, c, sigma2=1 )
{
    tau2 <- ifelse( lambda > sqrt( c )*sigma2,
                    2*(lambda + sigma2)^2 * (1 - c*(sigma2/lambda)^2),
                    0 )
    tau2
}

# Scaled limit of the standard deviation of a sample eigenvalue.  See
# "eigen.val.var".
#
eigen.val.sd <- function( lambda, c, sigma2=1 ) { 
    tau2 <- eigen.val.var( lambda, c, sigma2 )
    tau  <- sqrt( tau2 )
    tau
}

# Almost-sure limit of the square of the dot product between the sample 
# and the population eigenvector for a given eigenvalue.
#
# Arguments
# ---------
#     lambda : the population eigenvalue
#          c : the aspect ratio ( # of variables / sample size )
#     sigma2 : the noise variance
#
# Returns
# -------
#     r2 : the almost-sure limit of the square of the dot product
#
eigen.vec.cor2.mean <- function( lambda, c, sigma2=1 )
{
    r2 <- ifelse( lambda > sqrt( c )*sigma2,
                  ( lambda^2 - c*sigma2^2 )/( lambda*( lambda + c*sigma2 ) ),
                  0 )
    r2
}

# Almost-sure limit of the square of the length of the residual from the
# sample eigenvector projected onto the population eigenvector for a given
# eigenvalue.
#
# Arguments
# ---------
#     lambda : the population eigenvalue
#          c : the aspect ratio ( # of variables / sample size )
#     sigma2 : the noise variance
#
# Returns
# -------
#     rbar2 : the almost-sure limit of the square of length of the residual
#
eigen.vec.cor2.compl.mean <- function( lambda, c, sigma2=1 )
{
    rbar2 <- 1 - eigen.vec.cor2.mean( lambda, c, sigma2 )
    rbar2
}

# Shrink the sample eigenvalue to get a consistent estimate of the population
# eigenvalue.
#
# Arguments
# ---------
#          l: the sample eigenvalue
#          c: the aspect ratio ( # of variables / sample size )
#     sigma2: the noise variance
#
# Returns
# -------
#     lambda: a consistent estimate of the population eigenvalue when
#            $l > (1 + \sqrt(c))*\sigma^2$, and 'NA' otherwise.
#
shrink.eigen.val <- function( l, c, sigma2=1 )
{
    ow <- options( "warn" ); options( warn=-1 )
    lambda <- 0.5*( l - (1+c)*sigma2 
                  + sqrt( (l - (1+c)*sigma2 )^2 - 4*c*sigma2^2 ) )
    options( ow )
    
    ifelse( l > sigma2*(1 + sqrt( c ))^2, lambda, NA )
}

# Scaled limit of the variance of a shrunk sample eigenvalue.
#
# Arguments
# ---------
#     lambda: the population eigenvalue
#          c: the aspect ratio ( # of variables / sample size )
#     sigma2: the noise variance
#
# Returns
# -------
#     tau2.shrunk : the scaled variance of the shrunk sample eigenvalue.  
#                  This value divided by sample size is equal to the variance
#                  of the shrunk sample eigenvalue (with error o(1/n), where n
#                  is the sample size).
#
shrunk.eigen.val.var <- function( lambda, c, sigma2=1 )
{
    tau2.shrunk <- 2*( lambda + sigma2 )^2 / ( 1 - (c*sigma2/lambda)^2 ) 
    tau2.shrunk
}

# Scaled limit of the standard deviation of a shrunk sample eigenvalue.
# See "shrunk.eigen.val.var".
#
shrunk.eigen.val.sd <- function( lambda, c, sigma2=1 ) { 
    tau2.shrunk <- shrunk.eigen.val.var( lambda, c, sigma2 )
    tau.shrunk  <- sqrt( tau2.shrunk )
    tau.shrunk
}
