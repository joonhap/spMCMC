lplot <- function(...) plot(..., type='l')

algo <- 'spNUTS2'; eps <- 0.01; astar <- 0.65; repNo <- 1
load(file=paste('data/mvnorm_',algo,'_eps_',signif(eps,digits=2),'_astar_',astar,'_repNo_',repNo,'.RData',sep=''))

lplot(ESS.coda,ylim=c(0,max(ESS.coda)))
abline(h=c(0,Ma),lty=2)

par(mfrow=c(2,2))
plot(xxResult[,1])
plot(xxResult[,22])
plot(xxResult[,49])
plot(xxResult[,100])

cex <- 1.5
par(mfrow=c(1,2), mar=c(2,3,3,1),oma=c(0,2,0,0))
hist(unlist(result$log2len), main='', xlab='', ylab='', xaxt='n', yaxt='n', cex=cex)
mtext(text="Frequency", side=2, line=3, cex=cex)
mtext(text=bquote("log"[2]*"(trajectory length)"), side=3, line=0, cex=cex)
axis(side=1, at=min(unlist(result$log2len)):max(unlist(result$log2len)), cex.axis=cex); axis(side=2, cex.axis=cex)
hist(unlist(result$sym), main='', xlab='', ylab='', xaxt='n', yaxt='n', cex=cex)
mtext(text=bquote("Symmetry condition met"), side=3, line=0, cex=cex)
axis(side=1, at=c(0,1),cex.axis=cex); axis(side=2, cex.axis=cex)

par(mfrow=c(2,2))
plot(density(xxResult[,1])); abline(v=c(-sqrt(diag(Sigma)[1]),sqrt(diag(Sigma)[1])),col='red')
plot(density(xxResult[,10])); abline(v=c(-sqrt(diag(Sigma)[10]),sqrt(diag(Sigma)[10])),col='red')
plot(density(xxResult[,50])); abline(v=c(-sqrt(diag(Sigma)[50]),sqrt(diag(Sigma)[50])),col='red')
plot(density(xxResult[,100])); abline(v=c(-sqrt(diag(Sigma)[100]),sqrt(diag(Sigma)[100])),col='red')

par(mar=c(4,4,3,1)); cex=1.5
lplot(ESS.coda, ylim=c(0,max(ESS.coda)), xlab='', ylab='',xaxt='n',yaxt='n')
abline(h=c(0,M),lty=2)
axis(side=1,cex.axis=cex)
axis(side=2,cex.axis=cex)
mtext(text='component',side=1,line=2.5,cex=cex)
mtext(text='ESS',side=2,line=2.5,cex=cex)
mtext(text=bquote("min(ESS) = "*.(signif(minESS.coda,digits=3))*",  "*frac("min(ESS)","sec")*" = "*.(signif(minESS.coda.perSec,digits=3))), side=3, cex=cex)

par(mfrow=c(2,2), mar=c(2,2,2,1),oma=c(2,2,0,0))
acf(xxResult[,1], main='', xaxt='n', yaxt='n')
axis(side=1, cex.axis=cex); axis(side=2, cex.axis=cex)
mtext(text="1st component", side=3,line=0, cex=cex)
acf(xxResult[,12], main='', xaxt='n', yaxt='n')
axis(side=1, cex.axis=cex); axis(side=2, cex.axis=cex)
mtext(text="12th component", side=3,line=0, cex=cex)
acf(xxResult[,50], main='', xaxt='n', yaxt='n')
axis(side=1, cex.axis=cex); axis(side=2, cex.axis=cex)
mtext(text="50th component", side=3,line=0, cex=cex)
acf(xxResult[,100], main='', xaxt='n', yaxt='n')
axis(side=1, cex.axis=cex); axis(side=2, cex.axis=cex)
mtext(text="100th component", side=3,line=0, cex=cex)
mtext('lag',side=1,line=0.5,cex=cex, outer=TRUE)
mtext('ACF',side=2,line=0.5,cex=cex, outer=TRUE)

