## Linear sparse model
## Y = X beta + epsilon
## beta = delta * theta (entrywise product); delta: indicator (0 or 1), theta: effect size
## Y depends on only a small fraction of the variables (for which delta is one).
p <- 100
n.true.var <- 10 # number of true covariates
effect.size <- 2 # true covariates ~ Laplace(rate=1/effect.size)
nobs <- 100 # number of data points
noise.lvl <- .5 # the standard deviation of epsilon

seed <- 9487562 # random seed for model generation
set.seed(seed)

true.var <- sort(sample(x=1:p, size=10, replace=FALSE))
delta <- (1:p %in% true.var)
beta <- numeric(p)
beta[true.var] <- sample(x=c(-1,1), size=n.true.var, replace=TRUE) * rexp(n=n.true.var, rate=1/effect.size)

X <- matrix(rnorm(n=nobs*p), nrow=nobs, ncol=p)
epsilon <- rnorm(n=nobs, mean=0, sd=noise.lvl)
Y <- X %*% beta + epsilon


### diagnostic plot to see how well delta is identifiable
##coeff <- lm(Y~X)$coefficients[2:(p+1)]
##plot(beta, coeff)

