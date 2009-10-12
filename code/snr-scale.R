# snr-scale.R

require("RColorBrewer")
require("plyr")

MakeSNRScale <- function(levels.db=seq(-2,4), 
                         palette=brewer.pal(num.levels, "YlOrRd")[num.levels:1]) {
  num.levels <- length(levels.db)
  ramp       <- colorRamp(palette, space=c("rgb"))

  res <- list(num.levels=num.levels,
               min.level=levels.db[1],
               max.level=levels.db[num.levels],
                  levels=levels.db,
                 palette=palette,
                    ramp=ramp)
  class(res) <- c("SNRScale", class(res))
  res
}

ColorOnSNRScale <- function(snr.scale, db) {
  with(snr.scale, {
    u <- (db - min.level)/(max.level - min.level)
    aaply(u, 1, function(r) {
       do.call(rgb, c(as.list(ramp(r)), maxColorValue=255))
    }) 
  })     
}

plot.SNRScale <- function(x, ..., add=FALSE, pos="left", bty="n", cex=1.5) {
  if (!add) {
    plot(c(0,0), c(1,1), t="n", axes=FALSE, xlab="", ylab="")
  }
  num.levels <- x$num.levels
  labels     <- aaply(format(x$levels), 1, function(l) {
                    paste(l, "dB", sep="")})[num.levels:1]
  fill       <- x$palette[num.levels:1]
  legend(pos, legend=labels, fill=fill, cex=1.5, bty=bty, ...)
}

kDefaultSNRScale <- MakeSNRScale()

