### multiple proposal hybrid Hamiltonian Monte Carlo that samples both discrete and continuous variables
mphybridHMC <- function(X, disc.var, M, dmoves, rdirdist, ddirdist, leapfrog=lf, jsize, jumps, nu, pref=1) {
   ## X: existing history of chain. The last column is taken as the starting condition for the current iteration. If a vector, coerced into a matrix
    ## state variables are composed of both discrete and continuous variables
    ## disc.var: either the number of discrete variables, or the logical vector indicating whether a component is discrete or not. If just number (say 'n'), the first n variables are considered discrete variables. If a logical vector, the length should match the length of the state variable, and TRUE means a discrete variable.
    ## M: the max number of iterations
    ## veldist: the distribution from which the velocity vector is drawn
    ## dmoves: the number of move steps for updating discrete variables
    ## rdirdist: function that draws direction vectors for discrete variables
    ## ddirdist: evaluator of the probability of drawing a certain direction vector
    ## leapfrog: deterministic update map
    ## jsize, jumps: arguments to the leapfrog function for continuous variables
    ## nu: the distribution from which N and L are drawn
    dimst <- ifelse(is.vector(X), length(X), dim(X)[2]) # dimension of the state space
    dvar <- c() # index of discrete variables
    if (length(disc.var) == 1 && disc.var!=FALSE) {
        dvar <- 1:disc.var
    } else if (length(disc.var) == dimst) {
        dvar <- which(disc.var)
    } else {
        stop("disc.var must be either the number of discrete variables or a logical vector of the same length as the state variable, indicating whether each variable is discrete or not.")
    }
    cvar <- which(!(1:dimst %in% dvar)) # index of continuous variables
    muHMC <- function() {
        c(dmoves=dmoves, jsize=jsize, jumps=jumps)
    }
    rveldist <- function() { # uniform random directions for discrete variables and d-dimensional standard normal for continuous variables
        result <- numeric(dimst)
        result[dvar] <- rdirdist()
        result[cvar] <- rnorm(n=length(cvar))
        result
    }
    dveldist <- function(v, give_log=FALSE) {
        logcvardens <- sum(dnorm(v[cvar], log=TRUE))
        ifelse(give_log, log(ddirdist(v[dvar]))+logcvardens, ddirdist(v[dvar])*exp(logcvardens))
    }
    teo <- function(x, v, tau) {
        dmovRes <- dmove(x[dvar], v[dvar], tau["dmoves"])
        x[dvar] <- dmovRes[["x"]]
        v[dvar] <- dmovRes[["direc"]]
        cmovRes <- lf(x=x, v=v[cvar], jsize=tau["jsize"], jumps=tau["jumps"], fixed.var=dvar)
        x <- cmovRes[["x"]]
        v[cvar] <- cmovRes[["v"]]
        dmovRes <- dmove(x[dvar], v[dvar], tau["dmoves"])
        x[dvar] <- dmovRes[["x"]]
        v[dvar] <- dmovRes[["direc"]]
        list(x=x, v=v)
    }
    refop <- function(v) { -v }
    source('../mppd.R', local=TRUE) # load generic multiple proposal piecewise deterministic algorithm
    mppd(X, M, nu, muHMC, rveldist, dveldist, teo, refop, pref)
}



lf <- function(x, v, jsize, jumps, fixed.var=NULL) { # leapfrog algorithm
    ## jsize: jump size in each step
    ## jumps: number of jumps
    ## fixed.var: variables in x to be fixed in the leapfrog algorithm. Should be provided as a vector with elements in 1:length(x). The length of 'v' should equal the number of variables that are _not_ fixed.  The gradient function 'gd' should return a vector of the same length as 'v'.
    cvar <- which(!(1:length(x) %in% fixed.var)) # variables that are not fixed
    x[cvar] <- x[cvar] +.5*v*jsize
    j <- 1
    while (j < jumps) {
        v <- v + gd(x)*jsize
        x[cvar] <- x[cvar] + v*jsize
        j <- j+1
    }
    v <- v + gd(x)*jsize
    x[cvar] <- x[cvar] + .5*v*jsize

    list(x=x, v=v)
}


dmove <- function(x, direc, moves) { # move in discrete space
    ## direc: direction vector
    ## moves: number of moves
    for (n in 1:moves) {
        for (dd in 1:length(direc)) {
            if (direc[dd] > 0.5) {
                if (!is.max(x[dd], dd)) {
                    x[dd] <- x[dd]+1
                }
                else {
                    direc[dd] <- -direc[dd]
                }
            }
            else if (direc[dd] < -0.5) {
                if (!is.min(x[dd], dd)) {
                    x[dd] <- x[dd]-1
                }
                else {
                    direc[dd] <- -direc[dd]
                }
            }
        }
    }

    list(x=x, direc=direc)
}
