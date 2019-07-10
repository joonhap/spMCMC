## test the NUTS algorithm on the logistic regression
source('logitModel.R')
source('NUTS.R')

args = commandArgs(trailingOnly=TRUE)
#args <- c(0.002, .65, 1) # eps, astar, repNo

M <- 20200
burnout <- 200
Ma <- M-burnout
eps <- as.numeric(args[1]) # unstable if >~ .0027
algo <- 'NUTS'; adaptive.eps <- TRUE; adaptive.massInv <- TRUE
astar <- as.numeric(args[2])
repNo <- as.numeric(args[3])
set.seed(76234+578723*repNo)
print(paste(algo,'eps:',eps,'repNo',repNo))
x0 <- rnorm(25, sd=.1)
runtime <- system.time( {
    result <- NUTS(x0, M, eps, adaptive.eps=adaptive.eps, adaptive.massInv=adaptive.massInv, astar=astar)
}
)

## analysis of results
xxResult <- result$xx[1+burnout:M,]
nomove <- apply(xxResult[1+1:Ma,]-xxResult[1:Ma,], 1, function(x) all(x==0))
moveprob <- 1-sum(nomove)/Ma # the proportion that the MC makes a nonzero jump
print(paste('runtime:', runtime[1], 'moveprob', moveprob))
print(paste('offtrack', sum(is.na(unlist(result$log2len)))))

library(coda)
ESS.coda <- effectiveSize(xxResult)
minESS.coda <- min(ESS.coda)
avgESS.coda <- mean(ESS.coda)
minESS.coda.perSec <- minESS.coda/as.numeric(runtime)[1]*M/Ma
avgESS.coda.perSec <- avgESS.coda/as.numeric(runtime)[1]*M/Ma
print(paste('minESS', minESS.coda, 'minESS.perSec', minESS.coda.perSec, 'avgESS', avgESS.coda, 'avgESS.perSec', avgESS.coda.perSec))

log2len <- result$log2len; ontrack <- result$ontrack; eps.trace <- result$eps; accep.prob.trace <- result$accep.prob
save(eps, algo, repNo, runtime, moveprob, log2len, ontrack, eps.trace, accep.prob.trace, ESS.coda, minESS.coda, minESS.coda.perSec, avgESS.coda, avgESS.coda.perSec, file=paste('data/logit_',algo,'_eps_',signif(eps,digits=2),ifelse(adaptive.massInv,"","_fixedmass"),'_astar_',astar,'_repNo_',repNo,'.RData',sep=''))
