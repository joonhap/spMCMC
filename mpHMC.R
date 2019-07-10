### multiple proposal Hamiltonian Monte Carlo
mpHMC <- function(X, M, leapfrog=lf, jsize, jumps, nu, pref) {
    ## X: existing history of chain. The last column is taken as the starting condition for the current iteration. If a vector, coerced into a matrix
    ## M: the max number of iterations
    ## veldist: the distribution from which the velocity vector is drawn
    ## leapfrog: deterministic update map
    ## jsize, jumps: arguments to the leapfrog function
    ## nu: the distribution from which N and L are drawn
    muHMC <- function() {
        c(jsize=jsize*runif(n=1,min=0.8,max=1.2), jumps=jumps)
        ## vary leapfrog jump size randomly by +-20%
    }
    dimst <- ifelse(is.vector(X), length(X), dim(X)[1]) # dimension of the state space
    rveldist <- function() { # d dimensional standard normal
        rnorm(n=dimst)
    }
    dveldist <- function(v, give_log=FALSE) { 
        logresult <- sum(dnorm(v, log=TRUE))
        ifelse(give_log, logresult, exp(logresult))
    }
    teo <- function(x, v, tau) {
        lf(x, v, tau["jsize"], tau["jumps"])
    }
    refop <- function(x, v) { -v }
    source('../mppd.R', local=TRUE) # load generic multiple proposal piecewise deterministic algorithm
    mppd(X, M, nu, muHMC, rveldist, dveldist, teo, refop, pref=pref)
}



lf <- function(x, v, jsize, jumps) { # leapfrog algorithm
    ## jsize: the jump size in each step
    ## jumps: the number of jumps
    v <- v + .5*gd(x)*jsize
    j <- 1
    while (j < jumps) {
        x <- x + v*jsize
        v <- v + gd(x)*jsize
        j <- j+1
    }
    x <- x + v*jsize
    v <- v + .5*gd(x)*jsize

    list(x=x, v=v)
}

