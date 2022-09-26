#!/usr/bin/R
library(R.utils)
pdf("tmpexcursA_.pdf")
m <- numeric()
for( i in 1:countLines("./tmpexcurs")){
    x <- scan("./tmpexcurs", nlines=1, skip=i-1)
    m <- c( m, length(x) + 1)}
                                       plot(0,0, col='white', xlim=c(1, max(m)), ylim=c(0,1), las=1, xlab='time in generations', ylab='',  bty='l', cex.axis=1.4, cex.lab=1.5 )
#plot(0,0, col='white', xlim=c(0,1), ylim=c(0,1), las=1, xlab='time in generations', ylab='',  bty='l', cex.axis=1.4, cex.lab=1.5 )
for( i in 1:countLines("./tmpexcurs")){
    x <- c(  scan("./tmpexcurs", nlines=1, skip=i-1), 1)
 #                             lines(  (1:length(x))/length(x),   x, lwd=3, col=colors()[94+i], lty=i)}
   lines(   x, lwd=4, col=colors()[94+i], lty=i)}
graphics.off()
pdf("tmpexcursB_.pdf")
plot(0,0, col='white', xlim=c(0,1), ylim=c(0,1), las=1, xlab='', ylab='',  bty='l', cex.axis=1.4, cex.lab=1.5 )
for( i in 1:countLines("./tmpexcurs")){
    x <- c(  scan("./tmpexcurs", nlines=1, skip=i-1), 1)
                              lines(  (1:length(x))/length(x),   x, lwd=4, col=colors()[94+i], lty=i)}
print( c(  length(m), mean(m), sqrt(var(m)),  sqrt(var(m))/mean(m) ) )
graphics.off()
