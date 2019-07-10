## multiple proposal discretized bouncy particle sampler
mpdBPS <- function(X, M, ts, nu, refop1, refop2=NULL, mass=1, pref=0.2) {
    ## X: existing history of chain. The last column is taken as the initial condition. If a vector, coerced into a matrix
    ## refop1, refop2: reflection operators (refop2 is optional)
    ## mass: velocity distribution is defined as MVN(0,1/mass*identityMat)
    ## pref: the probability of refreshing the velocity vector
    if (is.vector(X)) { X <- cbind(X, deparse.level=0) } # make X a column vector without a label
    V <- cbind(rnorm(length(X[,1]),mean=0,sd=1/sqrt(mass)))
    for (i in 1:M) {
        NL <- nu()
        N <- NL[1] # the total number of proposals
        L <- NL[2] # the maximum number of accepted proposals
        Lam <- runif(1) # Uniform random number Lambda
        eps <- ts() # time step size
        Xc <- X[,dim(X)[2]] # current state of the chain
        Vc <- V[,dim(V)[2]] # current velocity vector
        Xn <- Xc # next state when proposals are rejected
        Vn <- refop1(Xc,Vc) # next velocity when proposals are rejected
        na <- 0 # number of accepted proposals
        Yc <- Xc # the sequence of proposals
        Wc <- Vc # the sequence of velocities
        baseline <- 1/tg(Xc) 
        for (n in 1:N) {
            Wn <- -refop1(Yc,Wc) # next velocity
            Yn <- Yc + Wn*eps # next proposal
            if (Lam < tg(Yn)*baseline) {
                na <- na + 1 }
            if (na == L) {
                Xn <- Yn
                Vn <- Wn
                break
            }
            Wc <- Wn
            Yc <- Yn
        }
        if(!is.null(refop2)) { Vn <- refop2(Xn,Vn) }
        if(runif(1) < pref) { Vn <- cbind(rnorm(length(X[,1]),mean=0,sd=1/sqrt(mass))) }
        X <- cbind(X, Xn, deparse.level=0)
        V <- cbind(V, Vn, deparse.level=0)
    }
    X
}
 

## reflection operators for multiple proposal discrete BPS
refop1 <- function(x,v) {
    gdU <- -gd(x) # gradient of U = -log pi(x)
    v - 2 * c(v%*%gdU) * gdU / c(gdU %*% gdU)
}
refop2 <- function(x,v) {
    gdU <- -gd(x) # gradient of U = -log pi(x)
    v - 2 * c(v%*%gdU) * gdU / c(gdU %*% gdU)
}
## time step selector for mpdBPS
ts <- function() {
    size <- .2
    rchisq(1, df=2) * size
}

