## test the mpNUTS algorithm on the 100 dimensional multivariate normal model
lplot <- function(...) plot(...,type='l')
source('mvnormModel.R')
source('HMC.R')
#source('mpHMC.R')
#source('NUTS.R')
#source('mpNUTS1.R')
#source('mpNUTS2.R')

library(mvtnorm)
library(coda)

epsVec <- c(0.008, 0.010, 0.012, 0.014)
reps <- 20
ess.coda.by.eps <- sapply(epsVec, function(eps) {
replicate(reps, expr={

eps <- 0.012 ## change this
x0 <- rmvnorm(n=1, sigma=Sigma)
M <- 1000
l <- 1 # a collection of l leapfrog jumps will be considered as a jumping unit
algo <- 'mpNUTS2'
runtime <- system.time(
    if (algo=='NUTS') {
        result <- NUTSresult <- NUTS(x0, M, eps)
    }
    else if (algo=='mpNUTS1') {
        N <- 5; jmax <- 15 # no more than 2^jmax jumps will be made from the kernel
        result <- mpNUTS1result <- mpNUTS1(x0, M, N, eps, l, jmax)
    }
    else if (algo=='mpNUTS2') {
        N <- 20; jmax <- 15 # no more than 2^jmax jumps will be made from the kernel
        result <- mpNUTS2result <- mpNUTS2(x0, M, N, eps, l, jmax)
    }
)

ESS.coda <- effectiveSize(result$xx)
minESS.coda <- min(ESS.coda)
minESS.coda.per.sec <- minESS.coda/as.numeric(runtime)[1]
c(minESS.coda, minESS.coda.per.sec)
}
)
}
, simplify="array")
dimnames(ess.coda.by.eps) <- list(c("minESS","minESS.per.sec"),1:reps,epsVec
)
##save(ess.coda.by.eps, file="data/mpNUTS2_ess_by_eps.RData")

algo <- 'mpNUTS2'; rf.eps <- 'ran'; rf.cos <- 'ran'
pdf(paste('figures/mvnormModel/',algo,'_ESS_by_eps_',rf.eps,'_eps','_',rf.cos,'_cos','.pdf',sep=''),width=8,height=4)
cex <- 1.2
par(mfrow=c(1,2), mar=c(2,4,1,0.5),oma=c(2,0,0,0))
boxplot(ess.coda.by.eps[1,,], xlab='', ylab='', xaxt='n', yaxt='n')
mtext(text="minimum ESS", side=2, line=2.5, cex=cex)
axis(side=1, at=1:4, labels=epsVec, cex.axis=cex); axis(side=2, cex.axis=cex)
boxplot(ess.coda.by.eps[2,,], xlab='', ylab='', xaxt='n', yaxt='n')
mtext(text="minimum ESS per second", side=2, line=2.5, cex=cex)
mtext(text="mean leapfrog jump size", side=1, line=.5, cex=cex, outer=TRUE)
axis(side=1, at=1:4, labels=epsVec, cex.axis=cex); axis(side=2, cex.axis=cex)
dev.off()



## analysis of results
xxResult <- result$xx
nomove <- apply(xxResult[1+1:M,]-xxResult[1:M,], 1, function(x) all(x==0))
moveprob <- 1-sum(nomove)/M # the proportion that the MC makes a nonzero jump
print(paste('runtime:', runtime[1], 'moveprob', moveprob))

par(mfrow=c(2,2))
plot(xxResult[,1])
plot(xxResult[,22])
plot(xxResult[,49])
plot(xxResult[,100])

par(mfrow=c(2,2))
plot(density(xxResult[,1])); abline(v=c(-sqrt(diag(Sigma)[1]),sqrt(diag(Sigma)[1])),col='red')
plot(density(xxResult[,10])); abline(v=c(-sqrt(diag(Sigma)[10]),sqrt(diag(Sigma)[10])),col='red')
plot(density(xxResult[,50])); abline(v=c(-sqrt(diag(Sigma)[50]),sqrt(diag(Sigma)[50])),col='red')
plot(density(xxResult[,100])); abline(v=c(-sqrt(diag(Sigma)[100]),sqrt(diag(Sigma)[100])),col='red')

plot(density(rnorm(M)))

library(coda)
ESS.coda <- effectiveSize(xxResult)
(minESS.coda <- min(ESS.coda))
(minESS.coda.perSec <- minESS.coda/as.numeric(runtime)[1])

rf.eps <- 'ran'; rf.cos <- 'ran'

pdf(paste('figures/mvnormModel/',algo,'_ESS_',rf.eps,'_eps','_',rf.cos,'_cos','.pdf',sep=''),width=5,height=5)

par(mar=c(4,4,3,1)); cex=1.5
lplot(ESS.coda, ylim=c(0,max(ESS.coda)), xlab='', ylab='',xaxt='n',yaxt='n')
abline(h=c(0,M),lty=2)
axis(side=1,cex.axis=cex)
axis(side=2,cex.axis=cex)
mtext(text='component',side=1,line=2.5,cex=cex)
mtext(text='ESS',side=2,line=2.5,cex=cex)
mtext(text=bquote("min(ESS) = "*.(signif(minESS.coda,digits=3))*",  "*frac("min(ESS)","sec")*" = "*.(signif(minESS.coda.perSec,digits=3))), side=3, cex=cex)

dev.off()

pdf(paste('figures/mvnormModel/',algo,'_ACF_',rf.eps,'_eps','_',rf.cos,'_cos','.pdf',sep=''),width=8,height=8)
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
dev.off()

if (algo=='NUTS') {
    pdf(paste('figures/mvnormModel/',algo,'_diagnostics_',rf.eps,'_eps','_',rf.cos,'_cos','.pdf',sep=''),width=4,height=4)

    par(mar=c(2,3,3,1),oma=c(0,2,0,0))
    hist(unlist(result$log2len), main='', xlab='', ylab='', xaxt='n', yaxt='n', cex=cex)
    mtext(text="Frequency", side=2, line=3, cex=cex)
    mtext(text=bquote("log"[2]*"(trajectory length)"), side=3, line=0, cex=cex)
    axis(side=1, at=min(unlist(result$log2len)):max(unlist(result$log2len)), cex.axis=cex); axis(side=2, cex.axis=cex)

    dev.off()
}
if (algo=='mpNUTS1') {
    pdf(paste('figures/mvnormModel/',algo,'_diagnostics_',rf.eps,'_eps','_',rf.cos,'_cos','.pdf',sep=''),width=12,height=4)

    cex <- 1.5
    par(mfrow=c(1,3), mar=c(2,3,3,1),oma=c(0,2,0,0))
    hist(result$nfail, main = '', ylab='', xlab='',xaxt='n', yaxt='n', cex=cex)
    mtext(text="Number of rejected proposals", side=3, line=0, cex=cex)
    mtext(text="Frequency", side=2, line=3, cex=cex)
    axis(side=1, cex.axis=cex); axis(side=2, cex.axis=cex)
    hist(unlist(result$log2len), main='', xlab='', ylab='', xaxt='n', yaxt='n', cex=cex)
    mtext(text=bquote("log"[2]*"(trajectory length)"), side=3, line=0, cex=cex)
    axis(side=1, at=min(unlist(result$log2len)):max(unlist(result$log2len)), cex.axis=cex); axis(side=2, cex.axis=cex)
    hist(unlist(result$sym), main='', xlab='', ylab='', xaxt='n', yaxt='n', cex=cex)
    mtext(text=bquote("Symmetry condition met"), side=3, line=0, cex=cex)
    axis(side=1, at=c(0,1),cex.axis=cex); axis(side=2, cex.axis=cex)

    dev.off()
}
if (algo=='mpNUTS2') {
    pdf(paste('figures/mvnormModel/',algo,'_diagnostics_',rf.eps,'_eps','_',rf.cos,'_cos','.pdf',sep=''),width=8,height=4)
    par(mfrow=c(1,2), mar=c(2,3,3,1),oma=c(0,2,0,0))
    hist(unlist(result$log2len), main='', xlab='', ylab='', xaxt='n', yaxt='n', cex=cex)
    mtext(text="Frequency", side=2, line=3, cex=cex)
    mtext(text=bquote("log"[2]*"(trajectory length)"), side=3, line=0, cex=cex)
    axis(side=1, at=min(unlist(result$log2len)):max(unlist(result$log2len)), cex.axis=cex); axis(side=2, cex.axis=cex)
    hist(unlist(result$sym), main='', xlab='', ylab='', xaxt='n', yaxt='n', cex=cex)
    mtext(text=bquote("Symmetry condition met"), side=3, line=0, cex=cex)
    axis(side=1, at=c(0,1),cex.axis=cex); axis(side=2, cex.axis=cex)
    dev.off()
}    

