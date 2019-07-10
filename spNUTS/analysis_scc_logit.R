logit.spNUTS1.result <- function(eps, N, astar, repNo, M, burnout) {
    filename <- paste('data/logit_','spNUTS1','_eps_',signif(eps,digits=2),'_N_',N,'_astar_',astar,'_repNo_',repNo,'.RData',sep='')
    load(filename)
    return(c(as.numeric(runtime)[1], moveprob, mean(nfail), mean(unlist(log2len)), mean(unlist(sym)), minESS.coda, minESS.coda.perSec, avgESS.coda, avgESS.coda.perSec))
}
logit.spNUTS1.results <- function(eps, N, astar, repNoRange, M, burnout) {
    output <- t(sapply(repNoRange, function(x) logit.spNUTS1.result(eps, N, astar, x, M, burnout)))
    colnames(output) <- c('runtime', 'moveprob', 'nfail', 'log2len', 'sym', 'minESS', 'minESSps', 'avgESS', 'avgESSps')
    output
}
logit.spNUTS2.result <- function(eps, astar, repNo, M, burnout) {
    filename <- paste('data/logit_','spNUTS2','_eps_',signif(eps,digits=2),'_astar_',astar,'_repNo_',repNo,'.RData',sep='')
    load(filename)
    return(c(as.numeric(runtime)[1], moveprob, mean(unlist(log2len)), mean(unlist(sym),na.rm=TRUE), minESS.coda, minESS.coda.perSec, avgESS.coda, avgESS.coda.perSec))
}
logit.spNUTS2.results <- function(eps, astar, repNoRange, M, burnout) {
    output <- t(sapply(repNoRange, function(x) logit.spNUTS2.result(eps, astar, x, M, burnout)))
    colnames(output) <- c('runtime', 'moveprob', 'log2len', 'sym', 'minESS', 'minESSps', 'avgESS', 'avgESSps')
    output
}
logit.NUTS.result <- function(eps, astar, repNo, M, burnout) {
    filename <- paste('data/logit_','NUTS','_eps_',signif(eps,digits=2),'_astar_',astar,'_repNo_',repNo,'.RData',sep='')
    load(filename)
    return(c(as.numeric(runtime)[1], moveprob, mean(unlist(log2len)), minESS.coda, minESS.coda.perSec, avgESS.coda, avgESS.coda.perSec))
}
logit.NUTS.results <- function(eps, astar, repNoRange, M, burnout) {
    output <- t(sapply(repNoRange, function(x) logit.NUTS.result(eps, astar, x, M, burnout)))
    colnames(output) <- c('runtime', 'moveprob', 'log2len', 'minESS', 'minESSps', 'avgESS', 'avgESSps')
    output
}


eps <- 0.002; repNoRange <- 1:10; M <- 20200; burnout <- 200

logit_spNUTS1_.95_1 <- logit.spNUTS1.results(eps, 1, .95, repNoRange, M, burnout)
logit_spNUTS1_.95_5 <- logit.spNUTS1.results(eps, 5, .95, repNoRange, M, burnout)
logit_spNUTS1_.85_1 <- logit.spNUTS1.results(eps, 1, .85, repNoRange, M, burnout)
logit_spNUTS1_.85_5 <- logit.spNUTS1.results(eps, 5, .85, repNoRange, M, burnout)
logit_spNUTS1_.75_1 <- logit.spNUTS1.results(eps, 1, .75, repNoRange, M, burnout)
logit_spNUTS1_.75_5 <- logit.spNUTS1.results(eps, 5, .75, repNoRange, M, burnout)
logit_spNUTS1_.65_1 <- logit.spNUTS1.results(eps, 1, .65, repNoRange, M, burnout)
logit_spNUTS1_.65_5 <- logit.spNUTS1.results(eps, 5, .65, repNoRange, M, burnout)
logit_spNUTS1_.55_1 <- logit.spNUTS1.results(eps, 1, .55, repNoRange, M, burnout)
logit_spNUTS1_.55_5 <- logit.spNUTS1.results(eps, 5, .55, repNoRange, M, burnout)
logit_spNUTS1_.45_1 <- logit.spNUTS1.results(eps, 1, .45, repNoRange, M, burnout)
logit_spNUTS1_.45_5 <- logit.spNUTS1.results(eps, 5, .45, repNoRange, M, burnout)
##logit_spNUTS1_.35_1 <- logit.spNUTS1.results(eps, 1, .35, repNoRange, M, burnout)
##logit_spNUTS1_.35_5 <- logit.spNUTS1.results(eps, 5, .35, repNoRange, M, burnout)
logit_spNUTS2_.95 <- logit.spNUTS2.results(eps, .95, repNoRange, M, burnout)
logit_spNUTS2_.85 <- logit.spNUTS2.results(eps, .85, repNoRange, M, burnout)
logit_spNUTS2_.75 <- logit.spNUTS2.results(eps, .75, repNoRange, M, burnout)
logit_spNUTS2_.65 <- logit.spNUTS2.results(eps, .65, repNoRange, M, burnout)
logit_spNUTS2_.55 <- logit.spNUTS2.results(eps, .55, repNoRange, M, burnout)
logit_spNUTS2_.45 <- logit.spNUTS2.results(eps, .45, repNoRange, M, burnout)
##logit_spNUTS2_.35 <- logit.spNUTS2.results(eps, .35, repNoRange, M, burnout)
logit_NUTS_.95 <- logit.NUTS.results(eps, .95, repNoRange, M, burnout)
logit_NUTS_.85 <- logit.NUTS.results(eps, .85, repNoRange, M, burnout)
logit_NUTS_.75 <- logit.NUTS.results(eps, .75, repNoRange, M, burnout)
logit_NUTS_.65 <- logit.NUTS.results(eps, .65, repNoRange, M, burnout)
logit_NUTS_.55 <- logit.NUTS.results(eps, .55, repNoRange, M, burnout)
logit_NUTS_.45 <- logit.NUTS.results(eps, .45, repNoRange, M, burnout)
##logit_NUTS_.35 <- logit.NUTS.results(eps, .35, repNoRange, M, burnout)

minESSps <- c(
logit_NUTS_.95[,"minESSps"],
logit_NUTS_.85[,"minESSps"],
logit_NUTS_.75[,"minESSps"],
logit_NUTS_.65[,"minESSps"],
logit_NUTS_.55[,"minESSps"],
logit_NUTS_.45[,"minESSps"],
##logit_NUTS_.35[,"minESSps"],
logit_spNUTS1_.95_1[,"minESSps"],
logit_spNUTS1_.85_1[,"minESSps"],
logit_spNUTS1_.75_1[,"minESSps"],
logit_spNUTS1_.65_1[,"minESSps"],
logit_spNUTS1_.55_1[,"minESSps"],
logit_spNUTS1_.45_1[,"minESSps"],
##logit_spNUTS1_.35_1[,"minESSps"],
logit_spNUTS1_.95_5[,"minESSps"],
logit_spNUTS1_.85_5[,"minESSps"],
logit_spNUTS1_.75_5[,"minESSps"],
logit_spNUTS1_.65_5[,"minESSps"],
logit_spNUTS1_.55_5[,"minESSps"],
logit_spNUTS1_.45_5[,"minESSps"],
##logit_spNUTS1_.35_5[,"minESSps"],
logit_spNUTS2_.95[,"minESSps"],
logit_spNUTS2_.85[,"minESSps"],
logit_spNUTS2_.75[,"minESSps"],
logit_spNUTS2_.65[,"minESSps"],
logit_spNUTS2_.55[,"minESSps"],
logit_spNUTS2_.45[,"minESSps"]
##logit_spNUTS2_.35[,"minESSps"]
)

avgESSps <- c(
logit_NUTS_.95[,"avgESSps"],
logit_NUTS_.85[,"avgESSps"],
logit_NUTS_.75[,"avgESSps"],
logit_NUTS_.65[,"avgESSps"],
logit_NUTS_.55[,"avgESSps"],
logit_NUTS_.45[,"avgESSps"],
##logit_NUTS_.35[,"avgESSps"],
logit_spNUTS1_.95_1[,"avgESSps"],
logit_spNUTS1_.85_1[,"avgESSps"],
logit_spNUTS1_.75_1[,"avgESSps"],
logit_spNUTS1_.65_1[,"avgESSps"],
logit_spNUTS1_.55_1[,"avgESSps"],
logit_spNUTS1_.45_1[,"avgESSps"],
##logit_spNUTS1_.35_1[,"avgESSps"],
logit_spNUTS1_.95_5[,"avgESSps"],
logit_spNUTS1_.85_5[,"avgESSps"],
logit_spNUTS1_.75_5[,"avgESSps"],
logit_spNUTS1_.65_5[,"avgESSps"],
logit_spNUTS1_.55_5[,"avgESSps"],
logit_spNUTS1_.45_5[,"avgESSps"],
##logit_spNUTS1_.35_5[,"avgESSps"],
logit_spNUTS2_.95[,"avgESSps"],
logit_spNUTS2_.85[,"avgESSps"],
logit_spNUTS2_.75[,"avgESSps"],
logit_spNUTS2_.65[,"avgESSps"],
logit_spNUTS2_.55[,"avgESSps"],
logit_spNUTS2_.45[,"avgESSps"]
##logit_spNUTS2_.35[,"avgESSps"]
)

runtime <- c(
logit_NUTS_.95[,"runtime"],
logit_NUTS_.85[,"runtime"],
logit_NUTS_.75[,"runtime"],
logit_NUTS_.65[,"runtime"],
logit_NUTS_.55[,"runtime"],
logit_NUTS_.45[,"runtime"],
##logit_NUTS_.35[,"runtime"],
logit_spNUTS1_.95_1[,"runtime"],
logit_spNUTS1_.85_1[,"runtime"],
logit_spNUTS1_.75_1[,"runtime"],
logit_spNUTS1_.65_1[,"runtime"],
logit_spNUTS1_.55_1[,"runtime"],
logit_spNUTS1_.45_1[,"runtime"],
##logit_spNUTS1_.35_1[,"runtime"],
logit_spNUTS1_.95_5[,"runtime"],
logit_spNUTS1_.85_5[,"runtime"],
logit_spNUTS1_.75_5[,"runtime"],
logit_spNUTS1_.65_5[,"runtime"],
logit_spNUTS1_.55_5[,"runtime"],
logit_spNUTS1_.45_5[,"runtime"],
##logit_spNUTS1_.35_5[,"runtime"],
logit_spNUTS2_.95[,"runtime"],
logit_spNUTS2_.85[,"runtime"],
logit_spNUTS2_.75[,"runtime"],
logit_spNUTS2_.65[,"runtime"],
logit_spNUTS2_.55[,"runtime"],
logit_spNUTS2_.45[,"runtime"]
##logit_spNUTS2_.35[,"runtime"]
)

log2len <- c( # final value of adaptively tuned leapfrog jump size
logit_NUTS_.95[,"log2len"],
logit_NUTS_.85[,"log2len"],
logit_NUTS_.75[,"log2len"],
logit_NUTS_.65[,"log2len"],
logit_NUTS_.55[,"log2len"],
logit_NUTS_.45[,"log2len"],
##logit_NUTS_.35[,"log2len"],
logit_spNUTS1_.95_1[,"log2len"],
logit_spNUTS1_.85_1[,"log2len"],
logit_spNUTS1_.75_1[,"log2len"],
logit_spNUTS1_.65_1[,"log2len"],
logit_spNUTS1_.55_1[,"log2len"],
logit_spNUTS1_.45_1[,"log2len"],
##logit_spNUTS1_.35_1[,"log2len"],
logit_spNUTS1_.95_5[,"log2len"],
logit_spNUTS1_.85_5[,"log2len"],
logit_spNUTS1_.75_5[,"log2len"],
logit_spNUTS1_.65_5[,"log2len"],
logit_spNUTS1_.55_5[,"log2len"],
logit_spNUTS1_.45_5[,"log2len"],
##logit_spNUTS1_.35_5[,"log2len"],
logit_spNUTS2_.95[,"log2len"],
logit_spNUTS2_.85[,"log2len"],
logit_spNUTS2_.75[,"log2len"],
logit_spNUTS2_.65[,"log2len"],
logit_spNUTS2_.55[,"log2len"],
logit_spNUTS2_.45[,"log2len"]
##logit_spNUTS2_.35[,"log2len"]
)

astarValues <- c(0.95,0.85,0.75,0.65,0.55,0.45)
algorithms <- c('NUTS','spNUTS1 (N=1)','spNUTS1 (N=5)','spNUTS2')

## min ESS per second vs. astar and algorithms
gen_pdf <- TRUE; cex <- 1.3
if(gen_pdf) { pdf(file='figures/logitModel/minESSperSec_astar_algo.pdf', width=10, height=3) }
par(mar=c(4,4.5,.3,.3))
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),minESSps, pch=3, xaxt='n',yaxt='n', xlab='',ylab='', cex=1, ylim=c(0,max(minESSps)))
axis(side=2, cex.axis=cex)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cex)
axis(side=1, line=2, at=(1+length(astarValues))/2+length(astarValues)*0:(length(algorithms)-1), labels=algorithms, lty=0, cex.axis=cex)
mtext(side=2, line=3, text='minimum ESS per second', cex=cex)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
if(gen_pdf) { dev.off() }

## avg ESS per second vs. astar and algorithms
gen_pdf <- TRUE; cex <- 1.3
if(gen_pdf) { pdf(file='figures/logitModel/avgESSperSec_astar_algo.pdf', width=10, height=3) }
par(mar=c(4,4.5,.3,.3))
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),avgESSps, pch=3, xaxt='n',yaxt='n', xlab='',ylab='', cex=1, ylim=c(0,max(avgESSps)))
axis(side=2, cex.axis=cex)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cex)
axis(side=1, line=2, at=(1+length(astarValues))/2+length(astarValues)*0:(length(algorithms)-1), labels=algorithms, lty=0, cex.axis=cex)
mtext(side=2, line=3, text='avgimum ESS per second', cex=cex)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
if(gen_pdf) { dev.off() }


## min and avg ESS per second vs. astar and algorithms
gen_pdf <- TRUE; cexlab <- 1.2; cexaxis <- 1.6; cexalgo <- 1.8; cexpt <- 8
if(gen_pdf) { pdf(file='figures/logitModel/ESSperSec_astar_algo.pdf', width=10, height=6.5) }
par(mar=c(4,4.5,0.3,.3)); par(mfrow=c(3,1), oma=c(2.8,0,2,0))
## minimum ESS per second
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),minESSps, pch='.', xaxt='n',yaxt='n', xlab='',ylab='', xaxs="i", cex=cexpt, xlim=c(1-.5,(length(algorithms)*length(astarValues))+.5), ylim=c(min(minESSps[-((length(astarValues)-1)*length(repNoRange)+repNoRange)]), max(minESSps)), log='y')
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='minimum ESS /sec', cex=cexlab)
axis(side=3, line=-1, outer=TRUE, at=(1+length(astarValues))/2+length(astarValues)*0:(length(algorithms)-1), labels=algorithms, lty=0, cex.axis=cexalgo)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
## average ESS per second
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),avgESSps, pch='.', xaxt='n',yaxt='n', xaxs="i", xlab='',ylab='', cex=cexpt, xlim=c(1-.5,(length(algorithms)*length(astarValues))+.5), ylim=c(min(avgESSps[-((length(astarValues)-1)*length(repNoRange)+repNoRange)]), max(avgESSps)), log='y')
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='average ESS /sec', cex=cexlab)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
## runtime in seconds
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)), runtime, pch='.', xaxt='n',yaxt='n', xlab='',ylab='', xaxs="i", cex=cexpt, xlim=c(1-.5,(length(algorithms)*length(astarValues))+.5), log='y')
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='runtime', cex=cexlab)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
mtext(side=1, line=5, text='Target acceptance probability', cex=cexlab)
if(gen_pdf) { dev.off() }













## producing plots for NESS 2019
minESSps <- c(
logit_NUTS_.95[,"minESSps"],
logit_NUTS_.85[,"minESSps"],
logit_NUTS_.75[,"minESSps"],
logit_NUTS_.65[,"minESSps"],
logit_NUTS_.55[,"minESSps"],
logit_NUTS_.45[,"minESSps"],
##logit_NUTS_.35[,"minESSps"],
logit_spNUTS1_.95_5[,"minESSps"],
logit_spNUTS1_.85_5[,"minESSps"],
logit_spNUTS1_.75_5[,"minESSps"],
logit_spNUTS1_.65_5[,"minESSps"],
logit_spNUTS1_.55_5[,"minESSps"],
logit_spNUTS1_.45_5[,"minESSps"]
##logit_spNUTS1_.35_5[,"minESSps"],
)
astarValues <- c(0.95,0.85,0.75,0.65,0.55,0.45)
algorithms <- c('NUTS', 'spNUTS')
## min ESS per second vs. astar and algorithms
gen_pdf <- TRUE; cexlab <- 1; cexaxis <- 1; cexalgo <- 1
if(gen_pdf) { pdf(file='figures/logitModel/minESSperSec_astar_NUTS_spNUTS1.pdf', width=6, height=3) }
par(mar=c(4,4.5,1.5,.3))
plot(rep(1:(length(algorithms)*length(astarValues)),each=length(repNoRange)),minESSps, pch='.', xaxt='n',yaxt='n', xlab='',ylab='', cex=3, ylim=c(0,max(minESSps)))
axis(side=2, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=NA, cex.axis=cexaxis)
axis(side=1, at=1:(length(algorithms)*length(astarValues)), labels=rep(astarValues,length(algorithms)), cex.axis=cexaxis, lty=0, line=-.2, las=3)
mtext(side=2, line=3, text='minimum ESS /sec', cex=cexlab)
axis(side=3, line=-1, at=(1+length(astarValues))/2+length(astarValues)*0:(length(algorithms)-1), labels=algorithms, lty=0, cex.axis=cexalgo)
abline(v=.5+length(astarValues)*1:(length(algorithms)-1))
mtext(side=1, line=3, text='Target acceptance probability', cex=cexlab)
if(gen_pdf) { dev.off() }
