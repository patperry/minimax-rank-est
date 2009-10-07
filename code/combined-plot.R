
read.sim <- function( i ) {
    t          <- 1:1000
    prefix     <- paste( "sim", i, "-", sep="" )
    r.gen      <- read.table( paste( prefix, "gen-rank.txt", sep="" ),
                              header=FALSE, row.names=NULL )
    err.gen    <- read.table( paste( prefix, "gen-err.txt", sep="" ),
                              header=FALSE, row.names=NULL )
    r.kn       <- read.table( paste( prefix, "kn-rank.txt", sep="" ),
                              header=FALSE, row.names=NULL )
    err.kn     <- read.table( paste( prefix, "kn-err.txt", sep="" ),
                              header=FALSE, row.names=NULL )
    r.minmax   <- read.table( paste( prefix, "minmax-rank.txt", sep="" ),
                              header=FALSE, row.names=NULL )
    err.minmax <- read.table( paste( prefix, "minmax-err.txt", sep="" ),
                              header=FALSE, row.names=NULL )

    r               <- data.frame( t, r.gen, r.minmax, r.kn )
    colnames( r )   <- c("t", "gen", "minmax", "kn")
    err             <- data.frame( t, err.gen, err.minmax, err.kn )
    colnames( err ) <- c("t", "gen", "minmax", "kn")
    list( t=t, r=r, err=err )
}

plot.sim( sim, left=TRUE ) {
    if( left ) { par( mar=c(1,4,7,1) + 0.1 )
    } else     { par( mar=c(1,1,7,4) + 0.1 ) }

    palette( brewer.pal( 7, "YlOrRd" ) )
    plot( c(1,T), c(-60,60), 
          t="n", 
       xlab="", 
       ylab="",
       axes=FALSE )

    if( left ) { mtext( "Simulation 1", line=4.5, cex=2 ) 
    } else     { mtext( "Simulation 2", line=4.5, cex=2 ) } 
   
    box()
    if( left ) { 
        axis( 2, labels=TRUE );
        mtext( "Active Signal", side=2, line=4.5, cex=1.25) 
        mtext( "Frequency (Degrees)", side=2, line=2.5, cex=0.9) 
    } else { 
        axis( 4, labels=TRUE ) 
    }
    axis( 3, labels=TRUE )


    for( tc in t.change ) {
        abline( v=tc, lty=2, col="lightgray" )
    }

    for( i in 1:length( snr.db ) ) {
        lines( t, freq.deg[i,], col=5 - snr.db[i], lwd=3 )
    }

    if( left ) { par( mar=c(1,4,1,1) + 0.1 )
    } else     { par( mar=c(1,1,1,4) + 0.1 ) }

    palette( brewer.pal( 9, "Set1" ) )
    plot( gen$r, 
        t="n", 
     xlab="", 
     ylab="",
     ylim=c(0,4),
     axes=FALSE )
    box()
    if( left ) { axis( 2, labels=TRUE );
                 mtext( "Rank", side=2, line=4.5, cex=1.25 )} else { axis( 4, labels=TRUE ) } 
 
    for( tc in t.change ) {
        abline( v=tc, lty=2, col="lightgray" )
    }

    lines( t, gen$r, col=2, lwd=2 )
    lines( t, res$r, col=3, lwd=1 )

    if( left ) { par( mar=c(6,4,1,1) + 0.1 )
    } else     { par( mar=c(6,1,1,4) + 0.1 ) }


    plot( t, r.gen, t='n', 
        xlab="", 
        ylab="",
        ylim=c(0,4),
        axes=FALSE )
    box()
    if( left ) {
        axis( 2, labels=TRUE );
        mtext( "Subspace Approximation Error", side=2, line=4.5, cex=1.25 )
        mtext( "Squared Frobenius Norm", side=2, line=2.5, cex=0.9 )
    } else {
        axis( 4, labels=TRUE ) 

    }
    axis( 1, labels=TRUE )
    
    for( tc in t.change ) {
        abline( v=tc, lty=2, col="lightgray" )
    }

    lines( t, r.gen, col=2 )
    lines( t, r.est, col=3 )

    mtext( "Time", side=1, line=4, cex=1.25 )

    par( mar=c(4,4,4,4) + 0.1 )
}

par( oma=c(0,0,0,0) )
layout( matrix( c(0,0,0,1:9), 3, 4 ), 
        widths=c(1.5,4,4,1.5),
       heights=c(6,5,6) )

