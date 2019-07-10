### multiple proposal HMC
source('leapfrog.R')

spHMC <- function(x0, M, N, L, eps, l, massInv=1, adaptive.eps=FALSE, adaptive.massInv=FALSE, astar=.65) {
    ## x0: initial state
    ## M: length of markov chain
    ## N: maximum number of trials per iteration
    ## L: the number of acceptable candidates to be found
    ## eps: the leapfrog jump size
    ## l: candidates are obtained at every l leapfrog jumps
    ## massInv: the velocity components are allowed have different scaling rates. Each velocity component is scaled by the inverse of the square root of the ``mass'' corresponding to that component. The argument massInv is supplied as a vector of length spdim, where spdim is the space dimension. The velocity distribution is given by N(0,diag(massInv)), and the leapfrog velocity update is given by v <- v + .5*massInv*gradient_of_log_target_density*time_increment.
    ## adaptive: if TRUE, the algorithm will implement the adaptive approach similar to that proposed in Haario et al. 2001, but the scaling is conducted componentwise. The initial covariance matrix is taken equal to massInv.
    ## astar: target acceptance probability when eps is adaptively tuned
    x <- x0
    xx <- matrix(NA, nrow=M+1, ncol=length(x0)) # markov chain history
    xx[1,] <- x0
    nfail <- rep(NA,M)
    if (!is.vector(massInv) || !length(massInv)%in%c(1,spdim)) { stop("massInv should be a vector of length either one or spdim.") }
    if(!adaptive.eps) { eps.trace <- NULL; accep.prob.trace <- NULL }
    if(!adaptive.massInv) { massInv.trace <- NULL}
    m0 <- 100 # after m0 the algorithm will use tuned parameters 
    if(adaptive.eps) { eps.trace <- rep(NA, M); accep.prob.trace <- rep(NA, M) }
    if(adaptive.massInv) { sumX <- x0; sumX2 <- x0^2; massInv.trace <- matrix(NA, nrow=M, ncol=length(x0)) } # to compute sample variances. massInv is set to the sample variances of the components of x for m>=m0.
    m <- 0
    propTotal <- 0 # total accumulated number of proposed candidates
    while (m < M) {
        y <- x; w.raw <- rnorm(spdim)
        Del <- -log(runif(1)) # maximum allowable energy increase
        E0 <- .5*sum(w.raw^2) - tg(y) # initial energy
        w <- w.raw*sqrt(massInv)
        n <- 0
        na <- 0
        while (n < N) {
            if(adaptive.eps && n==0) {
                lf.onestep <- leapfrog(y, w, eps, 1, massInv)
                y.onestep <- lf.onestep$x; w.onestep <- lf.onestep$v
                E.onestep <- -tg(y.onestep)+.5*sum(w.onestep^2/massInv)
                accep.prob <- exp(min(E0-E.onestep,0))  # take a sample acceptance probability (using the state obtained by the first leapfrog jump)
            }
            eps.r <- eps*runif(1,0.8,1.2)
            lf.new <- leapfrog(y, w, eps.r, l, massInv)
            yprime <- lf.new$x; wprime <- lf.new$v
            Eprime <- .5*sum(wprime^2/massInv)-tg(yprime)
            if (Eprime <= E0 + Del) { na <- na+1 }
            if (na==L) {
                x <- yprime; v <- wprime
                break
            }
            y <- yprime; w <- wprime
            n <- n+1
        }
        m <- m+1
        xx[m+1,] <- x
        nfail[m] <- n
        if(adaptive.eps) {
            eps.trace[m] <- eps
            accep.prob.trace[m] <- accep.prob
            if(m>=m0) { eps <- eps * exp((accep.prob-astar)*m^-.7) # adapt leapfrog step size: log(eps_{m+1}) = log(eps_m) + 1/m^{0.7} * (acceptance_probability - target_acceptance_probability)
            }
        }
        if(adaptive.massInv) {
            sumX <- sumX + x; sumX2 <- sumX2 + x^2
            massInv.trace[m,] <- massInv
            if(m>=m0) { massInv <- (sumX2 - sumX^2/(m+1))/m + rep(1e-6,spdim) }
        }
    }
    return(list(xx=xx, nfail=nfail, eps=eps.trace, accep.prob=accep.prob.trace, massInv=massInv.trace))
}

