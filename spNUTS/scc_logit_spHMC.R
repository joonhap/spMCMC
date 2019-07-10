## test the spHMC algorithm on the logistic regression
source('logitModel.R')
source('spHMC.R')

args = commandArgs(trailingOnly=TRUE)

spdim <- 25
M <- 1500
burnout <- 500
Ma <- M-burnout
eps <- as.numeric(args[1]) #0.001
algo <- 'spHMC'
L <- as.numeric(args[2]) # the number of acceptable candidates to find
N <- L*10 # the maximum number of trials before discarding the generated trajectory
l <- as.numeric(args[3]) # the unit number of leapfrog jumps
repNo <- as.numeric(args[4])
set.seed(18738+98424*repNo)
print(paste(algo,'eps:',eps,'N',N,'L',L,'l',l,'repNo',repNo))
x0 <- rnorm(spdim, sd=.1)
runtime <- system.time( {
    result <- spHMC(x0, M, N, L, eps, l)
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
minESS.coda.perSec <- minESS.coda/as.numeric(runtime)[1]*M/Ma
print(paste('minESS', minESS.coda, 'minESS.perSec', minESS.coda.perSec))


eps.trace <- result$eps; accep.prob.trace <- result$accep.prob; nfail <- result$nfail
save(eps, algo, L, N, l, repNo, runtime, moveprob, nfail, eps.trace, accep.prob.trace, ESS.coda, minESS.coda, minESS.coda.perSec, file=paste('data/logit_algo_',algo,'_eps_',signif(eps,digits=2),'_L_',L,'_N_',N,'_l_',l,'_repNo_',repNo,'.RData',sep=''))

