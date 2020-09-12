## test the spNUTS1 algorithm on the multivariate normal distribution
source('mvnormModel.R')
source('spNUTS1.R')

args = commandArgs(trailingOnly=TRUE)
##args = c(.01, 1, 1, 5, .75, 1, 1) # eps, ran.eps, ran.cos, N, astar, adaptive.massInv, repNo

M <- 20200
burnout <- 200
Ma <- M-burnout
eps <- as.numeric(args[1])
ran.eps <- as.numeric(args[2])
ran.cos <- as.numeric(args[3])
l <- 1 # a collection of l leapfrog jumps will be considered as a jumping unit
algo <- 'spNUTS1'; adaptive.eps <- TRUE; adaptive.massInv <- as.numeric(args[6])
N <- as.numeric(args[4])
astar <- as.numeric(args[5])
repNo <- as.numeric(args[7])
set.seed(13587+72365*repNo)
print(paste(algo,'eps:',eps,'ran.eps',ran.eps,'ran.cos',ran.cos,'N',N,'astar',astar,'repNo',repNo))
x0 <- rnorm(spdim)*sqrt(diag(Sigma))
runtime <- system.time( {
    jmax <- 15 # no more than 2^jmax jumps will be made from the kernel
    result <- spNUTS1(x0, M, N, eps, l, jmax, ran.eps, ran.cos, adaptive.eps=adaptive.eps, adaptive.massInv=adaptive.massInv, astar=astar)
###    result <- spNUTS1(x0, M, N, eps, l, jmax, ran.eps, ran.cos, massInv=diag(Sigma))
}
)

## analysis of results
xxResult <- result$xx[1+burnout:M,]
nomove <- apply(xxResult[1+1:Ma,]-xxResult[1:Ma,], 1, function(x) all(x==0))
moveprob <- 1-sum(nomove)/Ma # the proportion that the MC makes a nonzero jump
print(paste('runtime:', runtime[1], 'moveprob', moveprob))

library(coda)
ESS.coda <- effectiveSize(xxResult)
minESS.coda <- min(ESS.coda)
avgESS.coda <- mean(ESS.coda)
minESS.coda.perSec <- minESS.coda/as.numeric(runtime)[1]*M/Ma
avgESS.coda.perSec <- avgESS.coda/as.numeric(runtime)[1]*M/Ma

library(mcmcse)
ESS.mcmcse <- ess(xxResult)
minESS.mcmcse <- min(ESS.mcmcse)
avgESS.mcmcse <- mean(ESS.mcmcse)
minESS.mcmcse.perSec <- minESS.mcmcse/as.numeric(runtime)[1]*M/Ma
avgESS.mcmcse.perSec <- avgESS.mcmcse/as.numeric(runtime)[1]*M/Ma
ESS.mcmcse.tukey <- ess(xxResult, method='tukey')
minESS.mcmcse.tukey <- min(ESS.mcmcse.tukey)
avgESS.mcmcse.tukey <- mean(ESS.mcmcse.tukey)
minESS.mcmcse.tukey.perSec <- minESS.mcmcse.tukey/as.numeric(runtime)[1]*M/Ma
avgESS.mcmcse.tukey.perSec <- avgESS.mcmcse.tukey/as.numeric(runtime)[1]*M/Ma

print(paste('minESS.coda', minESS.coda, 'minESS.coda.perSec', minESS.coda.perSec, 'avgESS.coda', avgESS.coda, 'avgESS.coda.perSec', avgESS.coda.perSec, 'minESS.mcmcse', minESS.mcmcse, 'minESS.mcmcse.perSec', minESS.mcmcse.perSec, 'avgESS.mcmcse', avgESS.mcmcse, 'avgESS.mcmcse.perSec', avgESS.mcmcse.perSec, 'minESS.mcmcse.tukey', minESS.mcmcse.tukey, 'minESS.mcmcse.tukey.perSec', minESS.mcmcse.tukey.perSec, 'avgESS.mcmcse.tukey', avgESS.mcmcse.tukey, 'avgESS.mcmcse.tukey.perSec', avgESS.mcmcse.tukey.perSec))



nfail <- result$nfail; log2len <- result$log2len; sym <- result$sym; eps.trace <- result$eps; accep.prob.trace <- result$accep.prob
save(eps, algo, repNo, runtime, moveprob, nfail, log2len, sym, eps.trace, accep.prob.trace, ESS.coda, minESS.coda, minESS.coda.perSec, avgESS.coda, avgESS.coda.perSec, ESS.mcmcse, minESS.mcmcse, minESS.mcmcse.perSec, avgESS.mcmcse, avgESS.mcmcse.perSec, ESS.mcmcse.tukey, minESS.mcmcse.tukey, minESS.mcmcse.tukey.perSec, avgESS.mcmcse.tukey, avgESS.mcmcse.tukey.perSec, file=paste('data/2020Jan/mvnorm_',algo,'_eps_',signif(eps,digits=2),ifelse(adaptive.massInv,"","_fixedmass"),'_N_',N,'_astar_',astar,'_repNo_',repNo,'.RData',sep=''))




