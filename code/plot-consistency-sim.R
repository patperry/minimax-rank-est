# plot-consistency-sim.R 
# ----------------------
require("RColorBrewer")

source("doa-sim.R")
source("doa-sim-params.R")
source("cov-time-series.R")
source("rank-est-oracle.R")
source("rank-est-htest.R")
source("rank-est-minimax.R")

frequency         <- c(0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16, 32)
num.frequencies   <- length(frequency)
num.reps          <- 50
htest.rank.err1   <- matrix(NA, num.reps, num.frequencies)
htest.rank.err2   <- matrix(NA, num.reps, num.frequencies)
minimax.rank.err1 <- matrix(NA, num.reps, num.frequencies)
minimax.rank.err2 <- matrix(NA, num.reps, num.frequencies)

for (i in seq_len(num.reps)) {
    for (j in seq_len(num.frequencies)) {
        set.seed(i, "Mersenne-Twister")
        sim1     <- SampleDOASim(KavcicYang1(dim=9, freq.snapshot=frequency[j]))
        cov.est1 <- with(sim1, WindowedCovEst(snapshot, length.window=45*frequency[j], time))
        htest1   <- EstimateRankHTestTimeSeries(cov.est1)
        minimax1 <- EstimateRankMinimaxTimeSeries(cov.est1)
        
        htest.rank.err1[i,j]   <- mean(abs(sim1$par$cov.signal$rank - htest1$cov.ts$rank))
        minimax.rank.err1[i,j] <- mean(abs(sim1$par$cov.signal$rank - minimax1$cov.ts$rank))


        set.seed(i, "Mersenne-Twister")
        sim2     <- SampleDOASim(KavcicYang2(dim=9, freq.snapshot=frequency[j]))
        cov.est2 <- with(sim2, WindowedCovEst(snapshot, length.window=45*frequency[j], time))
        htest2   <- EstimateRankHTestTimeSeries(cov.est2)
        minimax2 <- EstimateRankMinimaxTimeSeries(cov.est2)
         
        htest.rank.err2[i,j]   <- mean(abs(sim2$par$cov.signal$rank - htest2$cov.ts$rank))
        minimax.rank.err2[i,j] <- mean(abs(sim2$par$cov.signal$rank - minimax2$cov.ts$rank))
    }
}

write.table(htest.rank.err1,   file="../results/htest-rank-err-1.txt",   row.names=FALSE, col.names=FALSE)
write.table(htest.rank.err2,   file="../results/htest-rank-err-2.txt",   row.names=FALSE, col.names=FALSE)
write.table(minimax.rank.err1, file="../results/minimax-rank-err-1.txt", row.names=FALSE, col.names=FALSE)
write.table(minimax.rank.err2, file="../results/minimax-rank-err-2.txt", row.names=FALSE, col.names=FALSE)
