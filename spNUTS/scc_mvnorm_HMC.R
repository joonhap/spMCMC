## test the HMC algorithm on the multivariate normal distribution
source('mvnormModel.R')
source('HMC.R')

args = commandArgs(trailingOnly=TRUE)
#args = c(.01, 50, .65, 1, 1) # eps, l, astar, adaptive.massInv, repNo

M <- 20200
burnout <- 200
Ma <- M-burnout
eps <- as.numeric(args[1]) 
algo <- 'HMC'; adaptive.eps <- TRUE; adaptive.massInv <- as.numeric(args[4])
l <- as.numeric(args[2]) # the number of leapfrog jumps
astar <- as.numeric(args[3])
repNo <- as.numeric(args[5])
set.seed(48486+145397*repNo)
print(paste(algo,'eps:',eps,'l',l,'astar',astar,'repNo',repNo))
x0 <- rnorm(spdim, sd=.1)
runtime <- system.time( {
    result <- HMC(x0, M, eps, l, adaptive.eps=adaptive.eps, adaptive.massInv=adaptive.massInv, astar=astar)
#    result <- HMC(x0, M, 1, l, massInv=diag(Sigma))
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
print(paste('minESS', minESS.coda, 'minESS.perSec', minESS.coda.perSec, 'avgESS', avgESS.coda, 'avgESS.perSec', avgESS.coda.perSec))

eps.trace <- result$eps; accep.prob.trace <- result$accep.prob
save(eps, algo, repNo, runtime, moveprob, eps.trace, accep.prob.trace, ESS.coda, minESS.coda, minESS.coda.perSec, avgESS.coda, avgESS.coda.perSec, file=paste('data/mvnorm_',algo,'_eps_',signif(eps,digits=2),ifelse(adaptive.massInv,"","_fixedmass"),'_astar_',astar,'_repNo_',repNo,'.RData',sep=''))
