## the NUTS sampler (Algorithm 3 in Hoffman, Gelman 2014)
source('leapfrog.R')

## the inner product between two vectors (with respect to a given metric)
innerProd <- function (v1,v2,massInv=1) { sum(v1*v2/massInv) }

## a function that keeps expanding a binary tree until a U-turn is observed
buildtree <- function(x, v, E0, Del, dir, j, eps, massInv=1) {
    DelMax <- 1000 # if too off-track, terminate the trajectory early
    if (j==0) {
        lf.new <- leapfrog(x, v, dir*eps, 1, massInv)
        xprime <- lf.new$x; vprime <- lf.new$v
        Eprime <- -tg(xprime)+.5*sum(vprime^2/massInv)
        nprime <- (Eprime <= E0 + Del)
        sprime <- (Eprime < E0 + Del + DelMax)
        oprime <- sprime # whether extending was ended due to the trajectory being off-track
        return(list(xminus=xprime,vminus=vprime,xplus=xprime,vplus=vprime,xprime=xprime,nprime=nprime,sprime=sprime,oprime=oprime))
    }
    else {
        bt.new <- buildtree(x, v, E0, Del, dir, j-1, eps, massInv)
        xminus <- bt.new$xminus; vminus <- bt.new$vminus; xplus <- bt.new$xplus; vplus <- bt.new$vplus
        xprime <- bt.new$xprime; nprime <- bt.new$nprime; sprime <- bt.new$sprime; oprime <- bt.new$oprime
        if (sprime) {
            if (dir == -1) {
                bt.new <- buildtree(xminus, vminus, E0, Del, dir, j-1, eps, massInv)
                xminus <- bt.new$xminus; vminus <- bt.new$vminus
            }
            else {
                bt.new <- buildtree(xplus, vplus, E0, Del, dir, j-1, eps, massInv)
                xplus <- bt.new$xplus; vplus <- bt.new$vplus
            }
            xdprime <- bt.new$xprime; ndprime <- bt.new$nprime
            sdprime <- bt.new$sprime; odprime <- bt.new$oprime
            if ((nprime+ndprime)==0 || runif(1) < ndprime/(nprime+ndprime)) {
                xprime <- xdprime
            }
            sprime <- sdprime * (innerProd(xplus-xminus,vminus,massInv) >=0) * (innerProd(xplus-xminus,vplus,massInv) >=0)
            oprime <- oprime * odprime
            nprime <- nprime + ndprime
        }
        return(list(xminus=xminus,vminus=vminus,xplus=xplus,vplus=vplus,xprime=xprime,nprime=nprime,sprime=sprime,oprime=oprime))
    }
}

## no-U-turn sampler by Hoffman and Gelman (2014)
NUTS <- function(x0, M, eps, massInv=1, adaptive.eps=FALSE, adaptive.massInv=FALSE, astar=.65) {
    ## massInv: the velocity components are allowed have different scaling rates. Each velocity component is scaled by the inverse of the square root of the ``mass'' corresponding to that component. The argument massInv is supplied as a vector of length spdim, where spdim is the space dimension. The velocity distribution is given by N(0,diag(massInv)), and the leapfrog velocity update is given by v <- v + .5*massInv*gradient_of_log_target_density*time_increment.
    ## adaptive: if TRUE, the algorithm will implement the adaptive approach similar to that proposed in Haario et al. 2001, but the scaling is conducted componentwise. The initial covariance matrix is taken equal to massInv.
    ## astar: target acceptance probability when eps is adaptively tuned
    x <- x0
    xx <- matrix(NA, nrow=M+1, ncol=length(x0)) # Markov chain history
    xx[1,] <- x0
    log2len <- rep(NA, M); ontrack <- rep(NA, M)
    if (!is.vector(massInv) || !length(massInv)%in%c(1,spdim)) { stop("massInv should be a vector of length either one or spdim.") }
    if(!adaptive.eps) { eps.trace <- NULL; accep.prob.trace <- NULL }
    if(!adaptive.massInv) { massInv.trace <- NULL }
    m0 <- 100 # after m0 the algorithm will use tuned parameters 
    if(adaptive.eps) { eps.trace <- rep(NA, M); accep.prob.trace <- rep(NA, M) }
    if(adaptive.massInv) { sumX <- x0; sumX2 <- x0^2; massInv.trace <- matrix(NA, nrow=M, ncol=length(x0)) } # to compute sample variances. massInv is set to the sample variances of the components of x for m>=m0. 
    m <- 0
    while (m < M) {
        v.raw <- rnorm(spdim) # refresh velocity, spdim should be defined in the parent env
        Del <- -log(runif(1)) # energy allowance
        E0 <- -tg(x) + .5*sum(v.raw^2) # energy of the initial state
        v <- v.raw*sqrt(massInv)
        if(adaptive.eps) {
            lf.onestep <- leapfrog(x, v, eps, 1, massInv)
            x.onestep <- lf.onestep$x; v.onestep <- lf.onestep$v
            E.onestep <- -tg(x.onestep)+.5*sum(v.onestep^2/massInv)
            accep.prob <- exp(min(E0-E.onestep,0))  # take a sample acceptance probability (using the state obtained by the first leapfrog jump)
        }
        xminus <- x; xplus <- x; vminus <- v; vplus <- v; j <- 0; n <- 1; s <- 1; o <- 1
        while (s==1) {
            dir.j <- sample(c(-1,1),1)
            if (dir.j==-1) {
                bt.new <- buildtree(xminus, vminus, E0, Del, dir.j, j, eps, massInv)
                xminus <- bt.new$xminus; vminus <- bt.new$vminus; xprime <- bt.new$xprime
            }
            else {
                bt.new <- buildtree(xplus, vplus, E0, Del, dir.j, j, eps, massInv)
                xplus <- bt.new$xplus; vplus <- bt.new$vplus; xprime <- bt.new$xprime
            }
            nprime <- bt.new$nprime; sprime <- bt.new$sprime; oprime <- bt.new$oprime
            if (sprime==1 && runif(1) < nprime/n) { x <- xprime }
            n <- n + nprime
            s <- sprime * (innerProd(xplus-xminus,vminus,massInv)>=0) * (innerProd(xplus-xminus,vplus,massInv)>=0)
            o <- o * oprime
            j <- j+1
        }
        m <- m+1
        xx[m+1,] <- x
        log2len[m] <- j; ontrack[m] <- o
        if(adaptive.eps) {
            eps.trace[m] <- eps
            accep.prob.trace[m] <- accep.prob
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
    return(list(xx=xx, log2len=log2len, ontrack=ontrack, eps=eps.trace, accep.prob=accep.prob.trace, massInv=massInv.trace))
}
