## test the spHMC algorithm on the multivariate normal distribution
source('mvnormModel.R')
source('spHMC.R')

args = commandArgs(trailingOnly=TRUE)
#args = c(.01, 1, 50, 0.65, 1, 1) # eps, L, l, astar, adaptive.massInv, repNo

M <- 200200
burnout <- 200
Ma <- M-burnout
eps <- as.numeric(args[1]) 
algo <- 'spHMC'; adaptive.eps <- TRUE; adaptive.massInv <- as.numeric(args[5])
L <- as.numeric(args[2]) # the number of acceptable candidates to find
N <- L*10 # the maximum number of trials before discarding the generated trajectory
l <- as.numeric(args[3]) # the unit number of leapfrog jumps
astar <- as.numeric(args[4])
repNo <- as.numeric(args[6])
set.seed(18738+98424*repNo)
print(paste(algo,'eps:',eps,'N',N,'L',L,'l',l,'astar',astar,'repNo',repNo))
x0 <- rnorm(spdim, sd=.1)
runtime <- system.time( {
    result <- spHMC(x0, M, N, L, eps, l, adaptive.eps=adaptive.eps, adaptive.massInv=adaptive.massInv, astar=astar)
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

eps.trace <- result$eps; accep.prob.trace <- result$accep.prob; nfail <- result$nfail
save(eps, algo, repNo, runtime, moveprob, nfail, eps.trace, accep.prob.trace, ESS.coda, minESS.coda, minESS.coda.perSec, avgESS.coda, avgESS.coda.perSec, ESS.mcmcse, minESS.mcmcse, minESS.mcmcse.perSec, avgESS.mcmcse, avgESS.mcmcse.perSec, ESS.mcmcse.tukey, minESS.mcmcse.tukey, minESS.mcmcse.tukey.perSec, avgESS.mcmcse.tukey, avgESS.mcmcse.tukey.perSec, file=paste('data/2020Jan/mvnorm_',algo,'_eps_',signif(eps,digits=2),ifelse(adaptive.massInv,"","_fixedmass"),'_l_',l,'_astar_',astar,'_M_',M,'_repNo_',repNo,'.RData',sep=''))

