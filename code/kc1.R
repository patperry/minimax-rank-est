# kc1.R
# -----

n     <- 9
r.max <- 4
beta  <- 0.9
#N     <- (r.max+1)*n # 1/( 1 - beta )
N <- 45
#N <- 5/(1 - beta )

min.signal <- 0.6309573

T <- 1000
t <- 1:T

t0 <-    1
t1 <-  150
t2 <-  300
t3 <-  400
t4 <-  500
t5 <-  700
t6 <-  800
t7 <- 1000

freq.deg <- matrix( NA, 6, T )
snr.db   <- rep( NA, 6 )
freq.deg[ 1,t0:t2 ] <- seq(  35, 25, len=t2-t0+1 ); snr.db[1] <-  0
freq.deg[ 2,t0:t4 ] <- rep(  -5,     len=t4-t0+1 ); snr.db[2] <-  3
freq.deg[ 3,t1:t3 ] <- rep(  50,     len=t3-t1+1 ); snr.db[3] <- -2
freq.deg[ 4,t1:t3 ] <- rep( -30,     len=t3-t1+1 ); snr.db[4] <- -2
freq.deg[ 5,t5:t7 ] <- seq(  20,  5, len=t7-t5+1 ); snr.db[5] <-  3
freq.deg[ 6,t6:t7 ] <- seq( -15, -5, len=t7-t6+1 ); snr.db[6] <-  0

t.change <- c(t0, t1, t2, t3, t4, t5, t6, t7)

