## test the mpNUTS algorithm on the Stochastic Volatility model
lplot <- function(...) plot(...,type='l')
source('SVmodel.R')
#source('NUTS.R')
#source('mpNUTS1.R')
source('mpNUTS2.R')

spdim <- T+1
x0 <- c(rexp(1,rate=0.01), rep(log(sd(ldyi)),T))
head(x0)
M <- 1000
eps <- 0.001
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
runtime

## analysis of results
xxResult <- result$xx
nomove <- apply(xxResult[1+1:M,]-xxResult[1:M,], 1, function(x) all(x==0))
moveprob <- 1-sum(nomove)/M # the proportion that the MC makes a nonzero jump
print(paste('runtime:', runtime[1], 'moveprob', moveprob))
sum(is.na(result$log2len))

par(mfrow=c(2,2))
plot(xxResult[,1])
plot(xxResult[,22])
plot(xxResult[,49])
plot(xxResult[,100])

library(coda)
ESS.coda <- effectiveSize(xxResult)
(minESS.coda <- min(ESS.coda))
(minESS.coda.perSec <- minESS.coda/as.numeric(runtime)[1])
lplot(ESS.coda)

