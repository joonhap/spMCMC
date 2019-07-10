## compute ESS from a correlated collection of sample points
effSize.hoffman.vec <- function(x, mean.x, var.x, ac.cutoff) {
    ## x: a vector of correlated sample points. 
    ## mean.x, var.x: the mean and the variance of the distribution for the sample points, or estimates thereof
    ## ac.cutoff: the maximum lag for which the estimate of the autocorrelation will be used to compute the effective sample size (see Appendix A of Hoffman and Gelman 2014)
    M <- length(x)
    ac <- c()
    for (s in 1:(M-1)) {
        ac.s <- sum((x[(s+1):M]-mean.x)*(x[1:(M-s)]-mean.x))/(M-s)/var.x
        if (abs(ac.s) < ac.cutoff) { break }
        ac <- c(ac,ac.s)
    }
    n.ac <- length(ac)
    return(M/(1+2*sum((1-1:n.ac/M)*ac)))
}

effSize.hoffman <- function(xx, mean.xx, var.xx, ac.cutoff) {
    ## xx: a vector or a matrix of correlated sample points. If a matrix, the ESS will be computed for each column of xx
    ## mean.xx, var.xx: the mean and the variance of the distribution for the sample points, or estimates thereof (if a matrix, these should be vectors of the length equal to the number of columns in xx)
    ## ac.cutoff: the maximum lag for which the estimate of the autocorrelation will be used to compute the effective sample size (see Appendix A of Hoffman and Gelman 2014)
    if (is.vector(xx)) { return(effSize.hoffman.vec(xx, mean.xx, var.xx, ac.cutoff)) }
    if (is.matrix(xx)) {
        dim.x <- dim(xx)[2]
        return(sapply(1:dim.x, function(d) { effSize.hoffman.vec(xx[,d], mean.xx[d], var.xx[d], ac.cutoff) }))
    }
}
