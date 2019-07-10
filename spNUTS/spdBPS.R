## multiple proposal discretized bouncy particle sampler
spdBPS <- function(X, M, ts, nu, refop1, refop2=NULL, mass=1, pref) {
    ## X: existing history of chain. The last column is taken as the initial condition. If a vector, coerced into a matrix
    ## M: the number of iterations
    ## ts: time step length for each jump
    ## nu: the distribution for N and L
    ## refop1, refop2: reflection operators (refop2 is optional and not currently used)
    ## mass: velocity distribution is defined as MVN(0,1/mass*identityMat)
    ## pref: the probability of refreshing the velocity vector
    muBPS <- function() {
        ts*runif(n=1,min=0.8,max=1.2)
        ## vary jump size randomly by +-20%
    }
    dimst <- ifelse(is.vector(X), length(X), dim(X)[1]) # dimension of the state space
    rveldist <- function() { # d dimensional standard normal
        rnorm(n=dimst, mean=0, sd=sqrt(1/mass))
    }
    dveldist <- function(v, give_log=FALSE) { 
        logresult <- sum(dnorm(x=v, mean=0, sd=sqrt(1/mass), log=TRUE))
        ifelse(give_log, logresult, exp(logresult))
    }
    teo <- function(x, v, tau) {
        w <- -refop1(x,v) # next velocity
        y <- x + w*tau # next proposal
        list(x=y, v=w)
    }
    source('sppd.R', local=TRUE) # load generic multiple proposal piecewise deterministic algorithm
    sppd(X, M, nu, muBPS, rveldist, dveldist, teo, refop, pref)
}

