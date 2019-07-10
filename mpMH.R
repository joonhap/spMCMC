## multiple proposal Metropolis Hastings
mpMH <- function(X, M, rprop, dprop, nu) {
    ## X: existing history of chain. The last column is taken as the starting condition for the current iteration. If a vector, coerced into a matrix
    ## rprop: draw proposal
    ## dprop: compute the density of the proposal kernel
    if (is.vector(X)) { X <- cbind(X, deparse.level=0) } # make X a column vector without a label 
    for (i in 1:M) {
        NL <- nu()
        N <- NL[1] # the total number of proposals
        L <- NL[2] # the maximum number of accepted proposals
        Lam <- runif(1) # Uniform random number Lambda
        Xc <- X[,dim(X)[2]] # current newest proposal
        Xn <- Xc # the next state when proposals are rejected
        na <- 0 # number of accepted proposals
        Y <- cbind(Xc, deparse.level=0) # the sequence of proposals
        lcumprob <- -log(tg(Xc)) # log cumulative probability ratios, defined as  prod_j=1^n q(Y_{j-1}|Y_j) / [ tg(Xc) prod_j=1^n q(Y_j|Y_{j-1}) ]
        for (n in 1:N) {
            Yc <- Y[,n] # current newest proposal
            Yn <- rprop(Yc) # next proposal
            Y <- cbind(Y, Yn, deparse.level=0)
            lcumprob <- lcumprob + dprop(Yc, Yn, log=TRUE) - dprop(Yn, Yc, log=TRUE)
            if (log(Lam) < lcumprob + log(tg(Yn))) {
                na <- na + 1 }
            if (na == L) {
                Xn <- Yn
                break
            }
        }
        X <- cbind(X, Xn, deparse.level=0)
    }
    X
}



## MALA proposal kernel
tsMALA <- .5 # time step
rMALA <- function(x) { # random draw from MALA proposal kernel
    y <- x + gd(x) * tsMALA + rnorm(length(x), mean=0, sd=sqrt(2*tsMALA))
    y
}
dMALA <- function(y,x,log=FALSE) { # density of MALA proposal kernel
    if (log) { ret <- sum(dnorm((y - x - gd(x) * tsMALA), mean=0, sd=sqrt(2*tsMALA), log=TRUE)) }
    else { ret <- prod(dnorm((y - x - gd(x) * tsMALA), mean=0, sd=sqrt(2*tsMALA))) }
    ret
}

## random walk proposal kernel
js <- 2 # jump size
rRW <- function(x) { # random draw from random walk proposal kernel
    y <- x + rnorm(length(x), mean=0, sd=js)
    y
}
dRW <- function(y,x,log=FALSE) { # density of random walk proposal kernel
    if (log) { ret <- sum(dnorm((y - x), mean=0, sd=js, log=TRUE)) }
    else { ret <- prod(dnorm((y - x), mean=0, sd=js)) }
    ret
}
