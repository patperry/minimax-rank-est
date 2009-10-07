
radar.defs <- (function() {
    source.pos <- c( 1/2, 4/5 )
    target.pos <- matrix( c( c( 1/5, 6/7 ), c( 9/10, 3/5 ), c( 1/2, 1/5 ) ),
                          ncol=2,
                         byrow=TRUE )
    n.target   <- nrow( target.pos )

    xlim       <- c( 0, 1 )
    ylim       <- c( 0, 1 )
    
    if( require( "RColorBrewer", quietly=TRUE ) ) {
        source.col <- "black"
        target.col <- brewer.pal( max( n.target, 3 ), "Set1" )[ 1:n.target ]
    } else {
        source.col <- 1
        target.col <- 1 + 1:n.target
    }
    
    source.pch <- 16
    source.cex <- 2
    source.lty <- 1
    source.lwd <- 1.5

    target.pch <- rep( source.pch, n.target )
    target.cex <- rep( source.cex, n.target )
    target.lty <- rep( source.lty, n.target )
    target.lwd <- rep( source.lwd, n.target )
    
    return( new.env() )
})()
    

setup <- function( xlim=get( "xlim", radar.defs ), 
                   ylim=get( "ylim", radar.defs ), ... ) {
    par( mar=c( 0, 0, 0, 0 ) )
    plot( x=xlim, 
          y=ylim, 
       xlim=xlim, 
       ylim=ylim, 
          t="n", 
       axes=FALSE, 
       xlab="",
       ylab="", ... )
    
    box()
}

source <- function( pos=get( "source.pos", radar.defs ), 
                    col=get( "source.col", radar.defs ), 
                    pch=get( "source.pch", radar.defs ),
                    cex=get( "source.cex", radar.defs ),
                    ... ) {
    points( x=pos[ 1 ], 
            y=pos[ 2 ], 
          col=col, 
          pch=pch,
          cex=cex, ... )
}

circle <- function( x, y, r, col, lty, lwd, ... ) {
    symbols( x=x,
             y=y,
       circles=r,
        inches=FALSE,
           add=TRUE,
            fg=col,
           lty=lty,
           lwd=lwd, ... )
}

source.wave <- function( t,
                         pos=get( "source.pos", radar.defs ), 
                         col=get( "source.col", radar.defs ),
                         lty=get( "source.lty", radar.defs ),
                         lwd=get( "source.lwd", radar.defs ),
                         ... ) {
    circle( pos[ 1 ], pos[ 2 ], t, col, lty, lwd, ... )
}


target.wave <- function( t, 
                         target.pos=get( "target.pos", radar.defs ), 
                         target.col=get( "target.col", radar.defs ),
                         target.lty=get( "target.lty", radar.defs ),
                         target.lwd=get( "target.lwd", radar.defs ),
                         source.pos=get( "source.pos", radar.defs ), 
                         source.col=get( "source.col", radar.defs ),
                         source.lty=get( "source.lty", radar.defs ),
                         source.lwd=get( "source.lwd", radar.defs ),
                         ... ) {

    hit <- t.hit( source.pos, target.pos )
    n.target <- length( hit )
    for( i in 1:n.target ) {
        if ( t >= hit[ i ] ) {
            circle( target.pos[ i, 1 ], target.pos[ i, 2 ], t - hit[ i ],
                    target.col[ i ], target.lty[ i ], target.lwd[ i ] )
        }
    }

}

targets <- function( pos=get( "target.pos", radar.defs ), 
                     col=get( "target.col", radar.defs ), 
                     pch=get( "target.pch", radar.defs ),
                     cex=get( "target.cex", radar.defs ),
                     ... ) {
    points( x=pos[ ,1 ], 
            y=pos[ ,2 ], 
          col=col, 
          pch=pch,
          cex=cex, ... )
}





t.hit <- function( source.pos=get( "source.pos", radar.defs ),
                   target.pos=get( "target.pos", radar.defs ) ) {
    usr  <- par( "usr" )
    pin  <- par( "pin" )
    inpx <- ( pin[ 1 ] )/( usr[ 2 ] - usr[ 1 ] )
    inpy <- ( pin[ 2 ] )/( usr[ 4 ] - usr[ 3 ] )
    
    apply( target.pos, 1, 
        function( t ) { 
            sqrt( sum( ( ( source.pos - t )*c( inpx, inpy ) )^2 ) )/inpx
        } )
}

plot.radar <- function( t ) {
    events <- t.hit()
    t.max  <- max( events )*3.75
    setup()
    source()
    targets()
    source.wave( t*t.max )
    target.wave( t*t.max )
}

n.frame <- 40
times   <- seq( 0, 1, length=n.frame )
for( i in 1:n.frame ) {
    dn <- "../movies/radar"
    fn <- paste( i, ".pdf", sep="" )
    fn <- ifelse( i < 10, paste( "0", fn, sep=""), fn)
    pdf( paste( dn, fn, sep="/"), width=6, height=4 )
    plot.radar( times[ i ] )
    dev.off()
}
