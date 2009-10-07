# kc2.R

n     <- 9
r.max <- 3
#N     <- (r.max+1)*n
N     <- 45

T <- 1000
t <- 1:T

freq.deg <- matrix( NA, 3, T )
snr.db  <- rep( NA, 3 )
freq.deg[1,] <- seq( 35,  20, len=T ); snr.db[1] <-  4
freq.deg[2,] <- seq( 20,  35, len=T ); snr.db[2] <-  4
freq.deg[3,] <- seq(  5, -10, len=T ); snr.db[3] <- -2

t.change <- c(500)
