## a basic leapfrog function
leapfrog <- function(x, v, jsize, jumps, massInv=1) { # leapfrog algorithm
    ## jsize: the jump size in each step
    ## jumps: the number of jumps
    ## massInv: the velocity components are allowed have different scaling rates. Each velocity component is scaled by the inverse of the square root of the ``mass'' corresponding to that component. The argument massInv is supplied as a vector of length spdim, where spdim is the space dimension. The velocity distribution is given by N(0,diag(massInv)), and the leapfrog velocity update is given by v <- v + .5*massInv*gradient_of_log_target_density*time_increment.
    v <- v + .5*massInv*gd(x)*jsize
    j <- 1
    while (j < jumps) {
        x <- x + v*jsize
        v <- v + massInv*gd(x)*jsize
        j <- j+1
    }
    x <- x + v*jsize
    v <- v + .5*massInv*gd(x)*jsize
    return(list(x=x, v=v))
}

