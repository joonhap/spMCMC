### multiple proposal No-U-Turn sampler (the second version)
## in this version of the spNUTS, the leapfrog trajectory is repeatedly extended such that the _number of acceptable proposals_ equals increasing powers of 2.
## note the difference from spNUTS in spNUTS.R where the trajectory length is given by l*2^j, regardless of the number of acceptable proposals along the trajectory.
source('leapfrog.R')

## the cosine of the angle between two vectors (with respect to a given metric)
cosAngle <- function (v1,v2,massInv=1) { sum(v1*v2/massInv)/sqrt(sum(v1*v1/massInv)*sum(v2*v2/massInv)) }

## proposal kernel that keeps extending the leapfrog trajectory until the number of accepted proposals are doubled
doublingKernel2 <- function(x, v, Del, eps, l, N, jmax, minCosAngle, massInv=1) {
    ## x,v: initial position and velocity
    ## Del: the energy allowance (for determining acceptability)
    ## eps: the leapfrog jump size
    ## l: l leapfrog jumps will make a unit move
    ## N: the maximum number of leapfrog jumps until an acceptable state is found. If no acceptable states are found in the first N jump, the kernel returns the initial state and terminates.
    ## jmax: the leapfrog trajectory is extended until the maximum number of 2^jmax acceptable states are found, if the stopping condition (observing a U-turn) has not been satisfied until then
    ## minCosAngle: kernel stops when the cosine of the angle formed by the displacement and velocity falls below minCosAngle
    E0 <- -tg(x) + .5*sum(v^2/massInv) # initial energy
    x0 <- x; v0 <- v
    xaccep <- x0; vaccep <- v0 # last accepted position and velocity
    ## find the first acceptable state
    i <- 0; xtry <- xaccep; vtry <- vaccep
    accep.prob <- NA # we will use the acceptance probability of the first proposal for tuning the leapfrog jump size
    while (i<N) {
        leap.new <- leapfrog(xtry, vtry, eps, l, massInv)
        xtry <- leap.new$x; vtry <- leap.new$v
        Etry <- -tg(xtry) + .5*sum(vtry*vtry/massInv)
        if(i==0) { accep.prob <- exp(min(E0-Etry,0)) } 
        if (Etry <= E0+Del) {
            xaccep <- xtry; vaccep <- vtry
            break
        }
        i <- i+1
    }
    if (i==N) { return(list(y=x0, w=v0, log2len=0, sym=NA, accep.prob=accep.prob)) }
    ## repeatedly extend the trajectory while doubling the number of acceptable states until a U-turn is observed
    j <- 1; # the number of acceptables in the trajectory will be 2^j
    while (j<=jmax) {
        k <- 1;
        ## store checkpoints for the symmetry condition
        xcheck <- matrix(NA, nrow=j, ncol=spdim); xcheck[1,] <- xaccep
        vcheck <- matrix(NA, nrow=j, ncol=spdim); vcheck[1,] <- vaccep
        checkpoint <- 2^(j-2); cpinterval <- j-2 # next checkpoint and the current interval between checkpoints
        while (k <= 2^(j-1)) {
            i <- 0; xtry <- xaccep; vtry <- vaccep
            while (i<N) {
                leap.new <- leapfrog(xtry, vtry, eps, l, massInv)
                xtry <- leap.new$x; vtry <- leap.new$v
                Etry <- -tg(xtry) + .5*sum(vtry*vtry/massInv)
                if (Etry <= E0+Del) {
                    xaccep <- xtry; vaccep <- vtry
                    break
                }
                i <- i+1
            }
            if (i==N) { return(list(y=x0, w=v0, log2len=j, sym=NA, accep.prob=accep.prob)) }
            if (k==checkpoint) {
                xcheck[j-cpinterval,] <- xaccep; vcheck[j-cpinterval,] <- vaccep
                cpinterval <- cpinterval-1
                checkpoint <- checkpoint + 2^cpinterval # set up the next checkpoint
            }
            k <- k+1
        }
        ## see if the trajectory made a U-turn
        if (cosAngle(xaccep-x0,v0,massInv) <=minCosAngle || cosAngle(xaccep-x0,vaccep,massInv) <=minCosAngle) {
            xfinal <- xaccep; vfinal <- vaccep
            break
        }
        j <- j+1
    }
    ## see if the symmetry condition holds
    sym <- 1
    if (j>jmax) { j <- jmax; xfinal <- xaccep; vfinal <- vaccep } # the case where no U-turn is observed till j=jmax (rare)
    for (i in 1:j) {
        if (cosAngle(xfinal-xcheck[i,],vfinal,massInv) <=minCosAngle || cosAngle(xfinal-xcheck[i,],vcheck[i,],massInv) <=minCosAngle) {
            sym <- 0; break
        }
    }
    if (sym) {
        return(list(y=xfinal, w=vfinal, log2len=j, sym=sym, accep.prob=accep.prob))
    }
    else {
        return(list(y=x0, w=v0, log2len=j, sym=sym, accep.prob=accep.prob))
    }
}

## multiple proposal No-U-Turn sampler (2nd version)
spNUTS2 <- function(x0, M, N, eps, l, jmax, ran.eps=FALSE, ran.cos=FALSE, massInv=1, adaptive.eps=FALSE, adaptive.massInv=FALSE, astar=.65) {
    ## ran.eps; logical value indicating whether the leapfrog step size is randomized
    ## ran.cos; logical value indicating whether the stopping cosine is randomized
    ## massInv: the velocity components are allowed have different scaling rates. Each velocity component is scaled by the inverse of the square root of the ``mass'' corresponding to that component. The argument massInv is supplied as a vector of length spdim, where spdim is the space dimension. The velocity distribution is given by N(0,diag(massInv)), and the leapfrog velocity update is given by v <- v + .5*massInv*gradient_of_log_target_density*time_increment.
    ## adaptive: if TRUE, the algorithm will implement the adaptive approach similar to that proposed in Haario et al. 2001, but the scaling is conducted componentwise. The initial covariance matrix is taken equal to massInv.
    ## astar: target acceptance probability when eps is adaptively tuned
    x <- x0
    xx <- matrix(NA,nrow=M+1,ncol=length(x0)) # markov chain history
    xx[1,] <- x0
    log2len <- rep(NA,M); sym <- rep(NA,M) # the number of extensions of trajectory and whether the symmetry condition was satisfied for each iteration
    if(!adaptive.eps) { eps.trace <- NULL; accep.prob.trace <- NULL }
    if(!adaptive.massInv) { massInv.trace <- NULL }
    m0 <- 100 # after m0 the algorithm will use tuned parameters    
    if(adaptive.eps) { eps.trace <- rep(NA, M); accep.prob.trace <- rep(NA, M) } 
    if(adaptive.massInv) { sumX <- x0; sumX2 <- x0^2;  massInv.trace <- matrix(NA, nrow=M, ncol=length(x0)) } # to compute sample variances. massInv is set to the sample variances of the components of x for m>=m0.
    m <- 0
    while (m < M) {
        y <- x; w.raw <- rnorm(spdim) # spdim should be in the parent env
        Del <- -log(runif(1)) # maximum allowable energy increase
        minCosAngle <- ifelse(ran.cos, runif(1,0,1), 0)
        w <- w.raw*sqrt(massInv)
        prop.new <- doublingKernel2(y, w, Del, eps*ifelse(ran.eps, runif(1,.8,1.2), 1), l, N, jmax, minCosAngle, massInv)
        x <- prop.new$y
        m <- m+1
        xx[m+1,] <- x
        log2len[m] <- prop.new$log2len
        sym[m] <- prop.new$sym
        if(adaptive.eps) {
            eps.trace[m] <- eps
            accep.prob.trace[m] <- accep.prob <- prop.new$accep.prob
            if(m>=m0) { eps <- eps * exp((accep.prob-astar)*m^-.7) # adapt leapfrog step size: log(eps_{m+1}) = log(eps_m) + 1/m^{0.7} * (acceptance_probability - target_acceptance_probability)
            }
        }
        if(adaptive.massInv) {
            sumX <- sumX + x; sumX2 <- sumX2 + x^2
            massInv.trace[m,] <- massInv
            if(m>=m0) {
                massInv <- (sumX2 - sumX^2/(m+1))/m + rep(1e-6,spdim)
            }
        }
    }
    return(list(xx=xx, log2len=log2len, sym=sym, eps=eps.trace, accep.prob=accep.prob.trace, massInv=massInv.trace))
}
