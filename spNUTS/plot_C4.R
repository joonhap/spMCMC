## 2d model with "four C's"
source('C4.R')

settings <- matrix(NA, 16, 4) # columns: inC, N, pref, rtype
settings[1,] <- c(1,1,0.1,3)
settings[2,] <- c(1,10,0.1,3)
settings[3,] <- c(1,20,0.1,1)
settings[4,] <- c(1,20,0.1,2)
settings[5,] <- c(1,20,0.1,3)
settings[6,] <- c(1,20,0.1,4)
settings[7,] <- c(1,20,0,1)
settings[8,] <- c(1,20,0,2)
settings[9,] <- c(1,20,0,3)
settings[10,] <- c(1,20,0,4)
settings[11,] <- c(0,1,0.1,1)
settings[12,] <- c(0,1,0.1,3)
settings[13,] <- c(0,10,0.1,1)
settings[14,] <- c(0,10,0.1,3)
settings[15,] <- c(0,20,0.1,1)
settings[16,] <- c(0,20,0.1,3)

L0 <- 1; jumpsize <- 0.1; M <- 120000; expNo <- 0

## plot Markov chain samples
library(raster)
tgval <- sapply( seq(-2,2,.02), function(x) { sapply(rev(seq(-2, 2,.02)), function(y) tg(c(x,y))) } )

dataDir <- 'data/C4model'; plotDir <- 'figures/C4model'
plot.this <- function(rowno) {
    inC <- settings[rowno,1]; N0 <- settings[rowno,2]; pref <- settings[rowno,3]
    Rtype <- c('gradRef', 'negGradRef', 'negId', 'mix')[settings[rowno,4]]
    filename <- paste('C4(', ifelse(inC, 'in', 'out'), ')_BPS_N',N0,'L',L0, 'pref',round(pref,digits=2), 'eps', round(jumpsize,digits=3), 'Rtype_', Rtype, '_M', M, '_', expNo, sep='')
    load(paste(dataDir,'/',filename,'.RData',sep=''))
    give_pdf <- TRUE
    if(give_pdf) { pdf(paste(plotDir,'/',filename,'.pdf',sep='')) }
    par(mar=c(0,0,0,0))
    tempMap <- raster(tgval^ifelse(inC,1,-1), xmn = -sidelen/2, xmx = sidelen/2, ymn = -sidelen/2, ymx=sidelen/2)
    plot(tempMap,axes = FALSE,legend=FALSE, xlim=c(-2,2),ylim=c(-2,2))
    points(x=XBPS[1,], y=XBPS[2,], pch='.')
    lines(x=XBPS[1,1+4*0:99], y=XBPS[2,1+4*0:99],col='red', lwd=1.5)
    if(give_pdf) { dev.off() }
}

for (i in 1:16)
    plot.this(i)




## convergence plots (distance between the sample mean and the true mean)
distsMAT <- matrix(NA, M, dim(settings)[1])
runtimes <- numeric(dim(settings)[1])
for (rowno in 1:16) {
    inC <- settings[rowno,1]; N0 <- settings[rowno,2]; pref <- settings[rowno,3]
    Rtype <- c('gradRef', 'negGradRef', 'negId', 'mix')[settings[rowno,4]]
    filename <- paste('C4(', ifelse(inC, 'in', 'out'), ')_BPS_N',N0,'L',L0, 'pref',round(pref,digits=2), 'eps', round(jumpsize,digits=3), 'Rtype_', Rtype, '_M', M, '_', expNo, sep='')
    load(paste(dataDir,'/',filename,'.RData',sep=''))
    distsMAT[,rowno] <- dists; runtimes[rowno] <- as.numeric(time.taken,units='secs')
}

## L2 error by runtime and various values of N (inC)
give_pdf <- FALSE
exps <- c(1,2,5) # experiments to be plotted
if(give_pdf) { pdf(paste(plotDir,'/L2error_by_runtime_and_N_inC_pref0.1_negID_M120000.pdf',sep=''), width=8, height=2.8) }
par(mar=c(3,3,.3,.3))
cexlab <- 1; cexaxis <- 1; ltys <- c(1:3); cols <- c(1,1,1)
xlimit <- c(0,max(runtimes[exps])); ylimit <- c(0,max(distsMAT[,exps]))
plot(NA, xaxt='n', yaxt='n', xlab='',ylab='', xlim=xlimit, ylim=ylimit)
for (i in 1:length(exps)) {
    lines((1:M)/M*runtimes[exps[i]], distsMAT[,exps[i]], lty=ltys[i], col=cols[i], lwd=2)
}
abline(h=0, lty=5)
axis(side=2, line=-.3, cex.axis=cexaxis, lty=0); axis(side=2, labels=NA)
axis(side=1, line=-.3, cex.axis=cexaxis, lty=0); axis(side=1, labels=NA)
mtext(side=2, line=1.8, text=expression('L'[2]*' error'), cex=cexlab)
mtext(side=1, line=1.8, text='runtime (sec)', cex=cexlab)
library(latex2exp)
legend('topright', legend=c('N=1','10','20'), lty=ltys, col=cols, lwd=2)
if(give_pdf) { dev.off() }

## L2 error by runtime and various values of N (outC)
give_pdf <- FALSE
exps <- c(12,14,16) # experiments to be plotted
if(give_pdf) { pdf(paste(plotDir,'/L2error_by_runtime_and_N_outC_pref0.1_negID_M120000.pdf',sep=''), width=8, height=2.8) }
par(mar=c(3,3,.3,.3))
cexlab <- 1; cexaxis <- 1; ltys <- c(1:3); cols <- c(1,1,1)
xlimit <- c(0,max(runtimes[exps])); ylimit <- c(0,max(distsMAT[,exps]))
plot(NA, xaxt='n', yaxt='n', xlab='',ylab='', xlim=xlimit, ylim=ylimit)
for (i in 1:length(exps)) {
    lines((1:M)/M*runtimes[exps[i]], distsMAT[,exps[i]], lty=ltys[i], col=cols[i], lwd=2)
}
abline(h=0, lty=5)
axis(side=2, line=-.3, cex.axis=cexaxis, lty=0); axis(side=2, labels=NA)
axis(side=1, line=-.3, cex.axis=cexaxis, lty=0); axis(side=1, labels=NA)
mtext(side=2, line=1.8, text=expression('L'[2]*' error'), cex=cexlab)
mtext(side=1, line=1.8, text='runtime (sec)', cex=cexlab)
library(latex2exp)
legend('topright', legend=c('N=1','10','20'), lty=ltys, col=cols, lwd=2)
if(give_pdf) { dev.off() }

## L2 error by runtime, p^ref=0 or 0.1, various velocity reflection operators (inC)
give_pdf <- TRUE
##plotDir <- '~/Research/spMCMC/report/StatComp/revision/figures/C4model'
exps <- c(7,10,5,3,6) # experiments to be plotted
if(give_pdf) { pdf(paste(plotDir,'/L2error_by_runtime_and_Rtype_inC_N20_pref0_0.1_M120000.pdf',sep=''), width=8, height=2) }
par(mar=c(3,3,.3,.3))
cexlab <- 1; cexaxis <- 1; ltys <- c(2,3,1,2,3); cols <- c(2,2,1,1,1)
xlimit <- c(0,max(runtimes[exps])); ylimit <- c(0,.5)#max(distsMAT[,exps])
plot(NA, xaxt='n', yaxt='n', yaxs='i', xlab='',ylab='', xlim=xlimit, ylim=ylimit)
for (i in 1:length(exps)) {
    lines((1:M)/M*runtimes[exps[i]], distsMAT[,exps[i]], lty=ltys[i], col=cols[i], lwd=1)
}
##abline(h=0, lty=5)
axis(side=2, line=-.3, cex.axis=cexaxis, lty=0); axis(side=2, labels=NA)
axis(side=1, line=-.3, cex.axis=cexaxis, lty=0); axis(side=1, labels=NA)
mtext(side=2, line=1.8, text=expression('L'[2]*' error'), cex=cexlab)
mtext(side=1, line=1.8, text='runtime (sec)', cex=cexlab)
library(latex2exp)
legend('topright', legend=c(expression('gradient,'*' p'^'ref'*'=0'),expression('mix,'*' p'^'ref'*'=0'),expression('-identity,'*' p'^'ref'*'=0.1'),expression('gradient,'*' p'^'ref'*'=0.1'),expression('mix,'*' p'^'ref'*'=0.1')), lty=ltys, col=cols, lwd=1, ncol=2)
if(give_pdf) { dev.off() }
