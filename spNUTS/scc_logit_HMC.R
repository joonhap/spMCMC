## test the HMC algorithm on the logistic regression
source('logitModel.R')
source('HMC.R')

args = commandArgs(trailingOnly=TRUE)

spdim <- 25
M <- 1500
burnout <- 500
Ma <- M-burnout
eps <- as.numeric(args[1]) #0.001
algo <- 'HMC'
l <- as.numeric(args[2]) # the number of leapfrog jumps
repNo <- as.numeric(args[3])
set.seed(48486+145397*repNo)
print(paste(algo,'eps:',eps,'l',l,'repNo',repNo))
x0 <- rnorm(25, sd=.1)
runtime <- system.time( {
    result <- HMC(x0, M, eps, l)
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

save(eps, algo, l, repNo, result, runtime, moveprob, ESS.coda, minESS.coda, minESS.coda.perSec, file=paste('data/logit_algo_',algo,'_eps_',signif(eps,digits=2),'_l_',l,'_repNo_',repNo,'.RData',sep=''))

