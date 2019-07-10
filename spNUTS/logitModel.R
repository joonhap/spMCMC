## logistic regression model for German credit data (Frank and Asuncion, 2010, see reference in Hoffman and Gelman, 2014)
germanCredit <- as.matrix(read.table('data/german.data-numeric'))
spdim <- 25
cd <- germanCredit[,-25] # customer data
y <- -2*germanCredit[,25]+3 # 1 if given credit, -1 otherwise
yx <- cbind(y, diag(y)%*%cd)
s2 <- 100 # sigma^2

## target density
tg <- function(ab, give_log=TRUE) {
    ans <- -sum(log(1+exp(-yx %*% ab))) - .5*sum(ab^2)/s2
    ifelse(give_log, ans, exp(ans))
}

## gradient of log target density
gd <- function(ab) {
    c((1/(exp(c(yx %*% ab)) + 1)) %*% yx) - ab/s2
}

