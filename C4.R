## two dimensional distribution on a square. The level sets of the density looks like four letter "C"s as below.
## | C | C |
## | C | C |
## The space is divided into four sub-squares of equal size, and the center of "C" is at the center of each sub-square.
## Each "C" is a circle with an opening. The opening may not to the right--i.e., it may not look like a "C" but a rotated "C".

sidelen <- 4 # the length of a side of the square-shaped space
radius <- .8 # the radius of a "C"
centers <- sidelen/4* matrix(c(1,1,-1,1,-1,-1,1,-1), ncol=2, byrow=TRUE) # the coordinates of four "C"s
slitsize <- pi/3 # the opening size in "C" in radian
slitangle <- c(-.48*pi, -0.06*pi, -.7*pi, .53*pi) # the angle of openings
Csp <- t(sapply(1:4, function(Cn) { c( centers[Cn,1] + radius*cos(slitangle[Cn]+slitsize/2), centers[Cn,2] + radius*sin(slitangle[Cn]+slitsize/2) ) } ) ) ## the coordinates of one of the two end points of the four "C"s, the point at which we start writing a "C" (Csp: C starting point)
Cep <- t(sapply(1:4, function(Cn) { c( centers[Cn,1] + radius*cos(slitangle[Cn]-slitsize/2), centers[Cn,2] + radius*sin(slitangle[Cn]-slitsize/2) ) } ) ) ## the coordinates of one of the two end points of the four "C"s, the point at which we end writing a "C" (Cep: C ending point)

anglediff <- function(angle1, angle2) {
    ## the angle between from angle2 to angle 1 between -pi and pi (counterclockwise)
    ad <- angle1 - angle2
    if (ad > pi) {
        return(ad - 2*pi)
    }
    else if (ad < -pi) {
        return(ad + 2*pi)
    }
    return(ad)
}

dist2Cn <- function(x, y, Cn){
    ## squared distance to the specified C (1,2,3,4 in counterclockwise order)
    ## x, y: coordinates, Cn: which circle (1,2,3,4)
    ## returns the distance and the coordinate of the closest point
    anglefromCen <- atan2(y-centers[Cn,2], x-centers[Cn,1]) # angle from the center of C
    ad <- anglediff(anglefromCen, slitangle[Cn])
    if (abs(ad) > slitsize/2) {
        ## if the angle is not in the opening angle
        cp <- c(centers[Cn,1]+radius*cos(anglefromCen), centers[Cn,2]+radius*sin(anglefromCen)) # the coordinate of the closest point
    }
    else if (ad > 0) {
        ## if the point is closer to the "start point" of the "C"
        cp <- Csp[Cn,]
    }
    else {
        ## if the point is closer to the "end point" of the "C"
        cp <- Cep[Cn,]
    }
    return(list(d2 = (x-cp[1])^2+(y-cp[2])^2, cp = cp))
}

dist2C <- function(x, y){
    ## squared distance to the closest C
    Cn <- ifelse(x>0, ifelse(y>0, 1, 4), ifelse(y>0, 2, 3)) # quadrant
    anglefromCen <- atan2(y-centers[Cn,2], x-centers[Cn,1]) # angle from the center of the C in the quadrant the point is in
    if (abs(anglediff(anglefromCen, slitangle[Cn])) > slitsize/2) {
        ## if the angle does not correspond to the opening, the distance to the C in the quadrant should be the minimum distance
        cp <- c(centers[Cn,1]+radius*cos(anglefromCen), centers[Cn,2]+radius*sin(anglefromCen)) # the coordinate of the closest point
    }
    else {
        ## otherwise, compute distance to all four C's and take minimum
        dists <- lapply(1:4, function(CC) dist2Cn(x, y, CC))
        cp <- dists[[which.min(sapply(1:4, function(n) dists[[n]][["d2"]]))]][["cp"]]
    }
    return(list(d2 = (x-cp[1])^2+(y-cp[2])^2, cp=cp))
}


charlen <- .12 # characteristic length at which the kernel decreases to exp(-1). (larger charlen means slower decrease)
depth <- 15 # how deep the values of the kernel changes
inC <- TRUE # should the prob density inside the letter C be higher than outside?
kern <- function(d2) {
    ## the kernel that converts a squared distance to a real value (the negative log target density)
    ## d2: squared distance
    return ( - ifelse(inC,1,-1) * depth* exp(-d2^2/charlen^4) )
}
dkern <- function(d2) {
    ## derivative of function kern with respect to d2
    return( ifelse(inC,1,-1) * kern(d2) * 2 * d2 / charlen^4 )
}


### model functions ###
bslope <- 100 # slope at which the energy barrier increases outside the square
## target density
tg <- function(x, give_log=FALSE) {
    ## x: two dimensional coordinate
    if (length(x) != 2) { stop("The length of x should be two.") }
    if (abs(x[1])<=sidelen/2 && abs(x[2])<=sidelen/2) {
        ans <- ifelse( give_log, -kern(dist2C(x[1],x[2])[["d2"]]), exp(-kern(dist2C(x[1],x[2])[["d2"]])) )
    }
    else {
        ans <- ifelse( give_log, -bslope*max(abs(x)), exp(-bslope*max(abs(x))) )
    }
    ans
}

## gradient of log target density
gd <- function(x) {
    ## x: two dimensional coordinate
    if (length(x) != 2) { stop("The length of x should be two.") }
    d2cp <- dist2C(x[1],x[2])
    d2 <- d2cp[["d2"]]
    if (abs(x[1])<=sidelen/2 && abs(x[2])<=sidelen/2) {
        ans <- dkern(d2) * 2 * (x - d2cp[["cp"]])
    }
    else {
        ans <- -bslope*c(abs(x[1])>abs(x[2]), abs(x[1])<abs(x[2]))*sign(x)
    }
    ans
}

## rescaled gradient of log target density
## in case the gradient is less than the machine precision, the reflection wrt the gradient can give out a NaN. To prevent this, a rescaled gradient that is not going to be zero is returned.
rgd <- function(x) {
    ## x: two dimensional coordinate
    if (length(x) != 2) { stop("The length of x should be two.") }
    d2cp <- dist2C(x[1],x[2])
    d2 <- d2cp[["d2"]]
    if (abs(x[1])<=sidelen/2 && abs(x[2])<=sidelen/2) {
        ans <- ifelse(inC,1,-1) * (x-d2cp[["cp"]])
    }
    else {
        ans <- -c(abs(x[1])>abs(x[2]), abs(x[1])<abs(x[2]))*sign(x)
    }
    ans
}


## plot the density for debugging purposes
##tgval <- sapply(seq(-2,2,.05), function(x) { sapply(seq(-2,2,.05), function(y) tg(c(x,y))) } )
##cols = matrix(hcl(h=scales::rescale(tgval, c(0, 200))), nrow(tgval))
##grid::grid.raster(cols)
