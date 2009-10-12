# main.R
# ------
require("RColorBrewer")

source("doa-sim.R")
source("doa-sim-params.R")
source("cov-time-series.R")
source("rank-est-htest.R")
source("rank-est-minimax.R")


set.seed(0, "Mersenne-Twister")
sim1     <- SampleDOASim(KavcicYang1(dim=9, freq.snapshot=1))
cov.est1 <- with(sim1, WindowedCovEst(snapshot, length.window=45, time))
htest1   <- EstimateRankHTestTimeSeries(cov.est1)
minimax1 <- EstimateRankMinimaxTimeSeries(cov.est1)

set.seed(0, "Mersenne-Twister")
sim2     <- SampleDOASim(KavcicYang2(dim=9, freq.snapshot=1))
cov.est2 <- with(sim2, WindowedCovEst(snapshot, length.window=45, time))
htest2   <- EstimateRankHTestTimeSeries(cov.est2)
minimax2 <- EstimateRankMinimaxTimeSeries(cov.est2)


palette(brewer.pal(9, "Set1"))

par(oma=c(0,0,0,0))
layout(matrix(c(0,0,0,1:9), 3, 4), widths=c(1.5,4,4,1.5), heights=c(6,5,6))

par(mar=c(1,4,7,1) + 0.1)
plot(sim1$par, ylim=c(-60, 60), xlab="", ylab="", axes=FALSE)
box()
axis(2, labels=TRUE)
mtext("Active Signal",       side=2, line=4.5, cex=1.25) 
mtext("Frequency (Degrees)", side=2, line=2.5, cex=0.9) 
axis(3, labels=TRUE)

par(mar=c(1,4,1,1) + 0.1)
plot(sim1$par$time, sim1$par$cov.signal$rank,
     t="n", xlab="", ylab="", ylim=c(0,4), axes=FALSE)
lines(sim1$par$landmarks)
box()
axis(2, labels=TRUE);
mtext("Rank", side=2, line=4.5, cex=1.25)
lines(sim1$par$time, sim1$par$cov.signal$rank, col=1, lwd=2)
lines(htest1$time,   htest1$cov.ts$rank,       col=2)
lines(minimax1$time, minimax1$cov.ts$rank,     col=3)


plot(c(0,1), c(0,1))


par(mar=c(1,1,7,4) + 0.1)
plot(sim2$par, ylim=c(-60, 60), xlab="", ylab="", axes=FALSE)
box()
axis(3, labels=TRUE)
axis(4, labels=TRUE)


par(mar=c(1,1,1,4) + 0.1) 

plot(sim2$par$time, sim2$par$cov.signal$rank,
     t="n", xlab="", ylab="", ylim=c(0,4), axes=FALSE)
lines(sim2$par$landmarks)
box()
axis(4, labels=TRUE);
lines(sim2$par$time, sim2$par$cov.signal$rank, col=1, lwd=2)
lines(htest2$time,   htest2$cov.ts$rank,       col=2)
lines(minimax2$time, minimax2$cov.ts$rank,     col=3)


plot(c(0,1), c(0,1))

par(mar=c(4,0,4,0) +0.1 )
plot(kDefaultSNRScale)

plot(c(0,0), c(1,1), t="n", axes=FALSE, xlab="", ylab="" )
legend( "left",
    legend=c("Oracle", "Hyp. Test", "Minimax"),
       fill=1:3, 
       bty="n",
       cex=1.5 )

plot(c(0,0), c(1,1), t="n", axes=FALSE, xlab="", ylab="" )
legend( "left",
    legend=c("Oracle", "Hyp. Test", "Minimax"),
       fill=1:3, 
       bty="n",
       cex=1.5 )

