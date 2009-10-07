# driver.R

pdf( "../plots/simulations.pdf", width=12, height=14 )

par( oma=c(0,0,0,0) )
layout( matrix( c(0,0,0,1:9), 3, 4 ), 
        widths=c(1.5,4,4,1.5),
       heights=c(6,5,6) )

left <- TRUE
source( "kc1.R" )
source( "analyze.R" )
source( "plot.R" )

write.table( gen$r, file="sim1-gen-rank.txt", row.names=FALSE, col.names=FALSE )
write.table( res$r, file="sim1-rank.txt", row.names=FALSE, col.names=FALSE )
write.table( r.gen, file="sim1-gen-err.txt", row.names=FALSE, col.names=FALSE )
write.table( r.est, file="sim1-err.txt", row.names=FALSE, col.names=FALSE )

left <- FALSE
source( "kc2.R" )
source( "analyze.R" )
source( "plot.R" )

write.table( gen$r, file="sim2-gen-rank.txt", row.names=FALSE, col.names=FALSE )
write.table( res$r, file="sim2-rank.txt", row.names=FALSE, col.names=FALSE )
write.table( r.gen, file="sim2-gen-err.txt", row.names=FALSE, col.names=FALSE )
write.table( r.est, file="sim2-err.txt", row.names=FALSE, col.names=FALSE )

palette( brewer.pal( 7, "YlOrRd" ) )

par( mar=c(4,0,4,0) +0.1 )
plot(c(0,0), c(1,1), t="n", axes=FALSE, xlab="", ylab="" )
legend( "left", 
    legend=c("  4dB", "  3dB", "  2dB", "  1dB", "  0dB", "-1dB", "-2dB"), 
       fill=1:7, 
       bty="n",
       cex=1.5 )

palette( brewer.pal( 9, "Set1" ) )
plot(c(0,0), c(1,1), t="n", axes=FALSE, xlab="", ylab="" )
legend( "left",
    legend=c("Truth", "Estimate"),
       fill=2:3, 
       bty="n",
       cex=1.5 )

plot(c(0,0), c(1,1), t="n", axes=FALSE, xlab="", ylab="" )
legend( "left",
    legend=c("True Rank", "Estimated Rank"),
       fill=2:3, 
       bty="n",
       cex=1.5 )

dev.off()
