
library( RMTstat )
library( RColorBrewer )

pdf( "../plots/tw-density.pdf" )

palette( brewer.pal( 9, "Set1" ) )
par( mar=c(4, 4, 2, 4) + 0.1 )

plot( function(x) dtw(x, beta=2), 
      xlim=c(-5,5), 
      xlab="x",
      ylab="Density", 
       col=2 )
curve( dtw(x, beta=1), add=TRUE, col=1 )

legend( "topright", 
        legend=expression( beta == 1, beta == 2 ),
           lty=1,
           col=c(1,2),
           bty='n',
         inset=0.02 )

axis( 3, labels=FALSE )
axis( 4, labels=FALSE )

dev.off()
