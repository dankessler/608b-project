source('sampler.R')

## generate a dataset

n <- 20
K <- 2

z <- matrix(0,nrow=n, ncol=K)
z[1:10,1] <- 1
z[11:20,2] <- 1

th <- matrix(c(0,10,50,100),nrow=2)

A <- matrix(rnorm(n^2, mean=z %*% th %*% t(z), sd = 1),nrow=n,ncol=n)

## try fitting sbm.gibbs to it

z.init <- z # initialize at the truth
th.init <- th # initialize at truth again
mu.init <- c(1/2,1/2)
n.iters <- 100
K <- 2
mu0 <- 10
sig0 <- 10
alpha <- c(1,1)

results <- sbm.gibbs(A,z.init,th.init,mu.init,n.iters,K,mu0,sig0,alpha)

results$z.p[2,,]
results$z.p[,1,]
results$th.p[,1,1]
