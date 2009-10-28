# plot-consistency-sim.R
# ----------------------

require("RColorBrewer")

tryCatch({
    frequency       <- 

    htest.rank.err1   <- read.table("../results/htest-rank-err-1.txt") 
    htest.rank.err2   <- read.table("../results/htest-rank-err-2.txt") 
    minimax.rank.err1 <- read.table("../results/minimax-rank-err-1.txt") 
    minimax.rank.err2 <- read.table("../results/minimax-rank-err-2.txt") 
}, error=function(e) {
    source("consistency-sim.R")    
})

palette(brewer.pal(9, "Set1"))

PlotSim <- function(htest.rank.err, minimax.rank.err,
                    frequency=c(0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16, 32),
                    freq.jitter=0.15,
                    ...) {
    num.frequencies <- length(frequency)
    num.reps        <- nrow(htest.rank.err)
    if (ncol(htest.rank.err) != num.frequencies
        || nrow(minimax.rank.err) != num.reps
        || ncol(minimax.rank.err) != num.frequencies) {
        stop("Arguments 'htest.rank.err', 'minimax.rank.err', and 'frequency'",
             " have inconsistent dimensions: ",
             toString(dim(htest.rank.err)), ", ",
             toString(dim(minimax.rank.err)), ", and ",
             length(frequency), ".")
    }

    htest.freq   <- frequency - freq.jitter
    htest.est    <- apply(htest.rank.err, 2, mean)
    htest.sd     <- apply(htest.rank.err, 2, sd)

    minimax.freq <- frequency + freq.jitter
    minimax.est  <- apply(minimax.rank.err, 2, mean)
    minimax.sd   <- apply(minimax.rank.err, 2, sd)

    xlim <- range(htest.freq, minimax.freq)
    ylim <- range(htest.est   - htest.sd,   htest.est   + htest.sd,
                  minimax.est - minimax.sd, minimax.est + minimax.sd)

    plot(xlim, ylim, t="n", ...)

    plot.vals <- function(freq, est, sd, pch=16, cex=1.2, lty=2, lwd=2, ...) {
        lines(freq, est, lty=lty, lwd=1, ...)
        segments(freq, est - sd, freq, est + sd, lwd=lwd, ...)
        segments(freq - freq.jitter, est - sd,
                 freq + freq.jitter, est - sd, lwd=lwd, ...)
        segments(freq - freq.jitter, est + sd,
                 freq + freq.jitter, est + sd, lwd=lwd, ...)

        points(freq, est, pch=pch, cex=cex, ...)
    }

    plot.vals(frequency, htest.est,   htest.sd,   col=2)
    plot.vals(frequency, minimax.est, minimax.sd, col=3)
}

par(oma=c(0,0,0,0))
layout(cbind(0, 1, 0, 2, 3), widths=c(1.5,4,0.4,4,1.5))

par(mar=c(6,4,1,1) + 0.1)
PlotSim(htest.rank.err1, minimax.rank.err1, 
        xlab="", ylab="", axes=FALSE)
box()
axis(1, labels=TRUE)
mtext("Sampling Frequency",      side=1, line=2.5, cex=1.25)
mtext("Snapshots per Unit Time", side=1, line=4.5, cex=0.9)
axis(2, labels=TRUE)
mtext("Rank Estimation Error",          side=2, line=4.5, cex=1.25)
mtext("Mean Absolute Difference", side=2, line=2.5, cex=0.90)
axis(3, labels=FALSE)
axis(4, labels=FALSE)

par(mar=c(6,1,1,4) + 0.1)
PlotSim(htest.rank.err2, minimax.rank.err2,
        xlab="", ylab="", axes=FALSE)
box() 
axis(1, labels=TRUE)
mtext("Sampling Frequency",      side=1, line=2.5, cex=1.25)
mtext("Snapshots per Unit Time", side=1, line=4.5, cex=0.9)
axis(2, labels=TRUE)
axis(3, labels=FALSE)
axis(4, labels=FALSE)

par(mar=c(4,0,4,0) +0.1 )
plot(c(0,0), c(1,1), t="n", axes=FALSE, xlab="", ylab="" )
legend( "left",
    legend=c("Hyp. Test", "Minimax"),
       fill=2:3, 
       bty="n",
       cex=1.5 )


