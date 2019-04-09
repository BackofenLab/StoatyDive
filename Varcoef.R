


mean_cov <- seq(1,100)
p <- .5 
n <- 11
varcoef <- (1 / sqrt(mean_cov * (1-p)))/sqrt(n-1)
plot(x=mean_cov, y=varcoef)

mean_cov <- 10
p <- seq(0.01,1.0, 0.01)
n <- 11
varcoef <- (1 / sqrt(mean_cov * (1-p)))/sqrt(n-1)
plot(x=p, y=varcoef)

mean_cov <- 10
p <- 0.5
n <- seq(1,100)
varcoef <- (1 / sqrt(mean_cov * (1-p)))/sqrt(n-1)
plot(x=n, y=varcoef)





