### Gaussian mixture model
nmix <- 5
pmix <- c(.1,.3,.1,.05,.45) # mixture probability
mixid <- 1:nmix # mixture id
means <- 0:(nmix-1) # the mean of each mixture component
sds <- c(1,.5,.7,2,1) # standard deviation of mixture components

## target density
tg <- function(x, give_log=FALSE) {
    ## x is composed of (delta, theta) where delta is discrete and theta is continuous
    ifelse(give_log, log(pmix[x[1]]) + dnorm(x[2], mean=means[x[1]], sd=sds[x[1]], log=TRUE), pmix[x[1]] * dnorm(x[2], mean=means[x[1]], sd=sds[x[1]], log=FALSE))
}

## gradient of the log target density for continuous variable
gd <- function(x) {
    -(x[2]-means[x[1]])/(sds[x[1]]*sds[x[1]])
}

## test for maximum or minimum for discrete variable
is.max <- function(dsc, dvarID) {
    dsc == mixid[nmix]
}

is.min <- function(dsc, dvarID) {
    dsc == mixid[1]
}

## distribution for the direction for which the discrete variable will be updated
rdirdist <- function() {
    sample(x=c(-1, 0, 1), size=1)
}

ddirdist <- function(dsc) {
    1/3
}
