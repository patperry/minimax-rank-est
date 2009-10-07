# analyze.R
# ---------

set.seed( 0, "Mersenne-Twister" )


source( "generate.R" )
source( "init.R" )
source( "fit.R" )
source( "select.R" )

sigma2.0 <- 1

gen <- generate( n, freq.deg, snr.db )
x   <- gen$x


select <- select.rank

res <- fit( x, N, sigma2.0, select )

l <- loss( gen$r, gen$w, gen$lambda, res$w )
r <- apply(l, 2, which.min) - 1
r.est  <- l[ res$r+1 + (n+1)*(t-1) ]
r.gen  <- l[ gen$r+1 + (n+1)*(t-1) ]
r.true <- l[     r+1 + (n+1)*(t-1) ]
