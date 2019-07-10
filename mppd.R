### multiple proposal piecewise deterministic MCMC
mppd <- function(X, M, nu, mu, rveldist, dveldist, teo, refop, pref=1, refop2=NULL, prefop2=NULL) {
    ## X: existing history of chain. The last column is taken as the starting condition for the current iteration. If a vector, coerced into a matrix
    ## M: the max number of iterations
    ## nu: the distribution from which N and L are drawn
    ## mu: the distribution from which the time step length is drawn
    ## rveldist: simulator of the velocity distribution
    ## dveldist: density of the velocity distribution
    ## teo: time evolution operator
    ## refop: velocity reflection operator
    ## pref: velocity refreshment probability
    ## optional velocity reflection operator at the end of each iteration
    ## the probability of optional velocity reflection
    if (is.vector(X)) { X <- cbind(X, deparse.level=0) } # make X a column vector without a label
    Xc <- X[,dim(X)[2]] # set initial state of the Markov chain
    V <- rveldist() # initialize the velocity vector
    for (i in 1:M) {
        NL <- nu()
        N <- NL[1] # the total number of proposals
        L <- NL[2] # the maximum number of accepted proposals
        tau <- mu() # time step length
        Y <- Xc # proposals
        W <- V # velocity vector
        logLam <- log(runif(1)) # log of a uniform random number Lambda
        na <- 0 # number of accepted proposals
        V <- refop(Y, V) # set the next state of velocity when proposals are rejected
        logdenom <- tg(Xc, give_log=TRUE) + dveldist(W, give_log=TRUE) # log(pi(X^{(i)})) + log(psi(V^{(i)}))
        for (n in 1:N) {
            YW <- teo(Y,W,tau) # next proposal
            Y <- YW[["x"]]
            W <- YW[["v"]]
            lognum <- tg(Y, give_log=TRUE) + dveldist(W, give_log=TRUE) # log(pi(Y_n)) + log(psi(W_n))
            if (logLam < lognum - logdenom) {
                na <- na + 1 }
            if (na == L) {
                Xc <- Y
                V <- W
                break
            }
        }
        X <- cbind(X, Xc, deparse.level=0)
        if (!is.null(refop2) && runif(1) < prefop2) {
            V <- refop2(Xc, V)
        }
        if (runif(1) < pref) {
            V <- rveldist()
        }
    }
    X
}
