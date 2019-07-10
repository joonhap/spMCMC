## model parameters
dim <- 1 # space dimension
mode_sep <- 2 # the distance of separation between nodes (see below)
sig <- 1 # the s.d. of mixture components
modes <- .5*mode_sep*cbind(diag(dim), -diag(dim)) # mode coordinates
ncomp <- dim(modes)[2] # number of mixture components

## target density
tg <- function(x, give_log=FALSE) {
    if (length(x) != dim) { stop("The length of x should match the space dimension.") }
    require(mvtnorm)
    result <- sum(apply(modes, 2, function(mu) dmvnorm(x, mean=mu, sigma=sig^2*diag(dim))))/ncomp # mixture density
    ifelse(give_log, log(result), result)
}

## marginal target density (for model investigation purposes)
mg <- function(x, margin, give_log=FALSE) {
    ## x: scalar
    ## margin: dimension for which marginal density will be computed
    result <- sum(sapply(modes[margin,], function(mu) dnorm(x, mean=mu, sd=sig)))/ncomp
    ifelse(give_log, log(result), result)
}    

## gradient of log target density
gd <- function(x) {
    if (length(x) != dim) { stop("The length of x should match the space dimension.") }
    require(mvtnorm)
    lpdf_comp <- apply(modes, 2, function(mu) dmvnorm(x, mean=mu, sigma=sig^2*diag(dim), log=TRUE)) # log pdf corresponding to each mixture component
    off_lpdf_comp <- lpdf_comp - max(lpdf_comp) # offsetted to prevent NaN
    disp <- outer(x, rep(1,ncomp)) - modes # displacements (x - mode)
    c((-disp/sig^2)%*%exp(off_lpdf_comp) / sum(exp(off_lpdf_comp))) # grad log mixture density
}
