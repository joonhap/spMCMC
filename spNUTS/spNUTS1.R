### multiple proposal No-U-Turn sampler
source('leapfrog.R')

## the cosine of the angle between two vectors (with respect to a given metric)
cosAngle <- function (v1,v2,massInv=1) { sum(v1*v2/massInv)/sqrt(sum(v1*v1/massInv)*sum(v2*v2/massInv)) }

## a function that builds a trajectory of length 2^j * l * eps
buildTrajectory <- function(x, v, eps, l, j, massInv=1) {
    ## x,v: initial state and velocity
    ## eps: leapfrog jump size
    ## l: the unit number of leapfrog jumps
    ## j: 2^j collections of l leapfrog jumps will be made
    ## this function returns the n_k-th state and velocity vectors for k in 0:j where n_k = sum_{k'=1}^{k} 2^{j-k'}
    ## these intermediate states and velocities are used to check the symmetry condition in the proposal kernel
    xx <- matrix(NA, nrow=j+1, ncol=spdim); vv <- matrix(NA, nrow=j+1, ncol=spdim)
    xx[1,] <- x; vv[1,] <- v
    k <- 1
    while (k<=j) {
        xv.new <- leapfrog(x, v, eps, 2^{j-k}*l, massInv)
        x <- xv.new$x; v <- xv.new$v
        k <- k + 1
        xx[k,] <- x; vv[k,] <- v
    }
    xv.new <- leapfrog(x, v, eps, l, massInv)
    xprime <- xv.new$x; vprime <- xv.new$v
    return(list(xprime=xprime,vprime=vprime,xx=xx,vv=vv))
}

## proposal kernel that keeps doubling the leapfrog trajectories until a U-turn
doublingKernel1 <- function(x, v, eps, l, jmax, minCosAngle, E0, DelStop, massInv=1) {
    ## jmax: the maximum of trajectory created will be 2^jmax * l
    ## minCosAngle: kernel stops when the cosine of the angle formed by the displacement and velocity falls below minCosAngle
    ## E0: initial energy
    ## DelStop: stop extending the trajectory if the energy increase exceeds DelStop
    x0 <- x; v0 <- v
    j <- 0
    bt.new <- buildTrajectory(x, v, eps, l, 0, massInv)
    xprime <- bt.new$xprime; vprime <- bt.new$vprime
    while (cosAngle(v0,xprime-x0,massInv) > minCosAngle && cosAngle(vprime,xprime-x0,massInv) > minCosAngle && j<jmax) {
#        if (.5*sum(vprime^2/massInv)-tg(xprime) > E0 + DelStop) {
#            return (list(y=x0, w=v0, log2len=NA, sym=NA))
#        }
        bt.new <- buildTrajectory(xprime, vprime, eps, l, j, massInv)
        xprime <- bt.new$xprime; vprime <- bt.new$vprime
        xx <- bt.new$xx; vv <- bt.new$vv
        j <- j + 1
    }
    k <- 1
    s <- 1 ## checking whether the symmetry condition is met
    Elast <- .5*sum(vprime^2/massInv) - tg(xprime) # the energy of the last state
    while (k<=j) {
        if (cosAngle(vv[k,],xprime-xx[k,],massInv) <= minCosAngle || cosAngle(vprime,xprime-xx[k,],massInv) <=minCosAngle) {
            s <- 0
            break
        }
#        if (.5*sum(vv[k,]^2/massInv)-tg(xx[k,]) > Elast + DelStop) {
#            s <- -1 # if stopped by being off-track (wrt the level set), set s=-1
#            break
#        }
        k <- k + 1
    }
    if (s==1) {
        return(list(y=xprime, w=vprime, log2len=j, sym=s))
    }
    else {
        return(list(y=x0, w=v0, log2len=j, sym=s))
    }
}

## multiple proposal No-U-Turn sampler
spNUTS1 <- function(x0, M, N, eps, l, jmax, ran.eps=FALSE, ran.cos=FALSE, massInv=1, adaptive.eps=FALSE, adaptive.massInv=FALSE, astar=.65) {
    ## x0: initial state
    ## M: length of markov chain
    ## for this version, L=1 (i.e., the first acceptable proposal is taken as the next state of the markov chain)
    ## also, velocity is refreshed at every new proposal
    ## ran.eps; logical value indicating whether the leapfrog step size is randomized
    ## ran.cos; logical value indicating whether the stopping cosine is randomized
    ## massInv: the velocity components are allowed have different scaling rates. Each velocity component is scaled by the inverse of the square root of the ``mass'' corresponding to that component. The argument massInv is supplied as a vector of length spdim, where spdim is the space dimension. The velocity distribution is given by N(0,diag(massInv)), and the leapfrog velocity update is given by v <- v + .5*massInv*gradient_of_log_target_density*time_increment.
    ## adaptive: if TRUE, the algorithm will implement the adaptive approach similar to that proposed in Haario et al. 2001, but the scaling is conducted componentwise. The initial covariance matrix is taken equal to massInv.
    ## astar: target acceptance probability when eps is adaptively tuned
    x <- x0 
    xx <- matrix(NA, nrow=M+1, ncol=length(x0)) # markov chain history
    xx[1,] <- x0
    DelStop <- 1000 # terminate kernel if the energy increase is greater than DelStop
    nfail <- rep(NA,M); log2len <- vector("list",M); sym <- vector("list",M) # number of failed proposals, the lengths of the trajectories for the proposals (in log_2), whether the symmetry condition was satisfied for each proposal
    if(!adaptive.eps) { eps.trace <- NULL; accep.prob.trace <- NULL }
    if(!adaptive.massInv) { massInv.trace <- NULL }
    m0 <- 100 # after m0 the algorithm will use tuned parameters
    if(adaptive.eps) { eps.trace <- rep(NA, M); accep.prob.trace <- rep(NA, M) }
    if(adaptive.massInv) { sumX <- x0; sumX2 <- x0^2; massInv.trace <- matrix(NA, nrow=M, ncol=length(x0)); } # to compute sample variances. massInv is set to the sample variances of the components of x for m>=m0. 
    m <- 0
    propTotal <- 0 # total accumulated number of proposed candidates
    while (m < M) {
        y <- x; w.raw <- rnorm(spdim)
        Del <- -log(runif(1)) # maximum allowable energy increase
        E0 <- .5*sum(w.raw^2) - tg(y) # initial energy (tg: log target density, should be available in the parent environment)
        w <- w.raw*sqrt(massInv)
        if(adaptive.eps) {
            lf.onestep <- leapfrog(y, w, eps, 1, massInv)
            y.onestep <- lf.onestep$x; w.onestep <- lf.onestep$v
            E.onestep <- -tg(y.onestep)+.5*sum(w.onestep^2/massInv)
            accep.prob <- exp(min(E0-E.onestep,0))  # take a sample acceptance probability (using the state obtained by the first leapfrog jump)
        }
        n <- 0
        jj <- c() # the lengths of simulated trajectories
        ss <- c() # whether the symmetry condition was satisfied at each proposal
        while (n < N) {
            minCosAngle <- ifelse(ran.cos, runif(1,0,1), 0)
            prop.new <- doublingKernel1(y, w, eps*ifelse(ran.eps,runif(1,0.8,1.2),1), l, jmax, minCosAngle, E0, DelStop, massInv); propTotal <- propTotal+1
            yprime <- prop.new$y; wprime <- prop.new$w
            Eprime <- .5*sum(wprime^2/massInv)-tg(yprime)
            jj <- c(jj,prop.new$log2len); ss <- c(ss,prop.new$sym)
            if (.5*sum(wprime^2/massInv)-tg(yprime) <= E0 + Del) { # (yprime, wprime) is acceptable
                x <- yprime; v <- wprime
                break
            }
            y <- yprime
            rmvnormVec <- rnorm(spdim); w <- rmvnormVec*sqrt(massInv)*sqrt(sum(wprime^2/massInv)/sum(rmvnormVec^2)) # refresh the velocity while keeping the magnitude
            n <- n+1
        }
        m <- m+1
        xx[m+1,] <- x
        nfail[m] <- n; log2len[[m]] <- jj; sym[[m]] <- ss
        if(adaptive.eps) {
            eps.trace[m] <- eps
            accep.prob.trace[m] <- accep.prob
            if(m>=m0) { eps <- eps * exp((accep.prob-astar)*m^-.7) # adapt leapfrog step size: log(eps_{m+1}) = log(eps_m) + 1/m^{0.7} * (acceptance_probability - target_acceptance_probability)
            }
        }            
        if(adaptive.massInv) {
            sumX <- sumX + x; sumX2 <- sumX2 + x^2;
            massInv.trace[m,] <- massInv
            if(m>=m0) { massInv <- (sumX2 - sumX^2/(m+1))/m + rep(1e-6,spdim) }
        }
    }
    return(list(xx=xx, nfail=nfail, log2len=log2len, sym=sym, eps=eps.trace, accep.prob=accep.prob.trace, massInv=massInv.trace))
}


