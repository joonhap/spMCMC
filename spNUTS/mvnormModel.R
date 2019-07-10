## model parameters
Sigma <- diag((0.01*(1:100))^2) # variance
spdim <- dim(Sigma)[1] # space dimension
invSigma <- solve(Sigma) # precision matrix
isDiag <- all(Sigma - diag(diag(Sigma)) == 0)

## target density
tg <- function(x, give_log=TRUE) {
    if (length(x) != spdim) { stop("The length of x should match the space dimension.") }
    if (isDiag) {
        ans <- sum(sapply(1:spdim, function(margin) -.5*x[margin]*x[margin]/Sigma[margin,margin]))
    }
    else {
        require(mvtnorm)
        ans <- dmvnorm(x, mean=rep(0,spdim), sigma=Sigma, log=TRUE)
    }
    ifelse(give_log, ans, exp(ans))
}

## marginal target density (for model investigation purposes)
mg <- function(x, margin, give_log=TRUE) {
    ## x: scalar
    ## margin: dimension for which marginal density will be computed
    dnorm(x, mean=mu, sd=sqrt(diag(Sigma)[margin]), log=give_log)
}    

## gradient of log target density
gd <- function(x) {
    if (length(x) != spdim) { stop("The length of x should match the space dimension.") }
    if (isDiag) {
        ans <- -x/diag(Sigma)
    }
    else {
        ans <- c(-invSigma %*% x)
    }
    ans
}

