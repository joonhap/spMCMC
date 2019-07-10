### MODELS ###

## 2d model with "four C's"
source('C4.R')
##source('../C4.R') # the original version

### ALGORITHMS ###

## multiple proposal discrete bouncy particle sampler
source("spdBPS.R")


## algorithmic parameters
args <- commandArgs(trailingOnly=TRUE)
#args <- c(N=1, L=1, pref=0.1, js=0.1, rtype=3, M=30000, expno=0) # N, L, pref, jumpsize, M, expNo
N0 <- as.numeric(args[1]); L0 <- as.numeric(args[2])
nu0 <- function() { # draws N and L (randomly or deterministically)
    c('N'=N0, 'L'=L0)
}
pref <- as.numeric(args[3])

jumpsize=as.numeric(args[4]) # discrete jump size

Rtype <- c('gradRef', 'negGradRef', 'negId', 'mix')[as.numeric(args[5])] # the type of reflection operator (gradRef: reflection with respect to the hyperplane perpendicular to the gradient of negative log target density; negGradRef: (-1) times gradRef; negId: multiply (-1); mix: mix gradRef and -Id with equal probability)
refop <- function(x, v) {
    if (Rtype != 'negId') { rgdU <- rgd(x) } # (rescaled) gradient of negative log target density (such that it does not become less than the machine precision and become zero)
    switch(Rtype,
           gradRef = v-2*sum(v*rgdU)*rgdU/sum(rgdU*rgdU),
           negGradRef = -v+2*sum(v*rgdU)*rgdU/sum(rgdU*rgdU),
           negId = -v,
           mix = { if(runif(1)<.5) { v-2*sum(v*rgdU)*rgdU/sum(rgdU*rgdU) } else { -v } }
           )
}

## initialization
init <- c(-.15, .01) # an arbitrary value

## run MCMC
M <- as.numeric(args[6])

expNo <- as.numeric(args[7])
seed <- 718542L + 987234L*expNo + 878234L*L0 + 7847261L*N0 + round(8723481*pref) + round(7987938*jumpsize) + 143*M
set.seed(seed)

start.time <- Sys.time()
XBPS <- spdBPS(X=init, M=M-1, ts=jumpsize, nu=nu0, refop1=refop, pref=pref)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


## convergence of the mean to the true mean (the origin)
cumsumXBPS <- apply(XBPS, 1, cumsum)
meanXBPS <- cumsumXBPS/1:M
dists <- apply(meanXBPS, 1, function(x) sqrt(sum(x^2)))


## save results
dataDir <- 'data/C4model'
filename <- paste('C4(', ifelse(inC, 'in', 'out'), ')_BPS_N',N0,'L',L0, 'pref',round(pref,digits=2), 'eps', round(jumpsize,digits=3), 'Rtype_', Rtype, '_M', M, '_', expNo, sep='')
info <- paste('N=',N0,', L=',L0, ', M=', M,
              '\n2 dimensional model with four "C"s ', ifelse(inC, '(inside)', '(outside)'),
              '\nsquare side length ', sidelen,
              '\nC slit size ', round(slitsize/pi,digits=2), 'pi', '\nslit opening angles (in counterclockwise order)', slitangle[1]/pi, ' ', slitangle[2]/pi, ' ', slitangle[3]/pi, ' ', slitangle[4]/pi ,' pi',
              '\nkernel depth ', depth, ' kernel characteristic length ', charlen,
              '\ndiscrete jump size ', jumpsize, '+-20%',
              '\np_ref: ', pref,
              '\nreflection operator: ', Rtype,
              '\nexperiment number: ', expNo,
              '\nseed= ', seed,
              '\n', sep='')
save_results <- FALSE
if(save_results) { save(XBPS, dists, time.taken, info, file=paste(dataDir,'/',filename,'.RData',sep='')) }


## make a plot 
library(raster)
tgval <- sapply( seq(-2,2,.02), function(x) { sapply(rev(seq(-2, 2,.02)), function(y) tg(c(x,y))) } )
par(mar=c(0,0,0,0))
tempMap <- raster(tgval, xmn = -sidelen/2, xmx = sidelen/2, ymn = -sidelen/2, ymx=sidelen/2)
plot(tempMap,axes = FALSE,legend=FALSE, xlim=c(-2,2),ylim=c(-2,2))
points(x=XBPS[1,], y=XBPS[2,], pch='.')
lines(x=XBPS[1,1+4*0:99], y=XBPS[2,1+4*0:99],col='red', lwd=1.5)
