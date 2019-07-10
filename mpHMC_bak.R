### multiple proposal Hamiltonian Monte Carlo
mpHMC <- function(X, M, leapfrog=lf, jsize, jumps, nu) {
    ## X: existing history of chain. The last column is taken as the starting condition for the current iteration. If a vector, coerced into a matrix
    ## M: the max number of iterations
    ## veldist: the distribution from which the velocity vector is drawn
    ## leapfrog: deterministic update map
    ## jsize, jumps: arguments to the leapfrog function
    ## nu: the distribution from which N and L are drawn
    require(mvtnorm)
    if (is.vector(X)) { X <- cbind(X, deparse.level=0) } # make X a column vector without a label 
    for (i in 1:M) {
        NL <- nu()
        N <- NL[1] # the total number of proposals
        L <- NL[2] # the maximum number of accepted proposals
        logLam <- log(runif(1)) # log of a uniform random number Lambda
        Xc <- X[,dim(X)[2]] # current newest proposal
        W <- c(rmvnorm(n=1,mean=rep(0,length(Xc)))) # the velocity vector
        Xn <- Xc # the next state when proposals are rejected
        na <- 0 # number of accepted proposals
        Y <- Xc # proposals
        niHam <- log(tg(Xc)) - .5*W*W # negative initial Hamitonian
        for (n in 1:N) {
            YW <- leapfrog(Y, W, jsize, jumps) # next proposal
            Y <- YW[["x"]]
            W <- YW[["v"]]
            nHam <- log(tg(Y)) - .5*W*W # negative Hamiltonian
            if (logLam < nHam - niHam) {
                na <- na + 1 }
            if (na == L) {
                Xn <- Y
                break
            }
        }
        X <- cbind(X, Xn, deparse.level=0)
    }
    X
}



lf <- function(x, v, jsize, jumps) {
    ## jsize: the jump size in each step
    ## jumps: the number of jumps
    x <- x+.5*v*jsize
    j <- 1
    while (j < jumps) {
        v <- v + gd(x)*jsize
        x <- x + v*jsize
        j <- j+1
    }
    v <- v + gd(x)*jsize
    x <- x + .5*v*jsize

    list(x=x, v=v)
}


