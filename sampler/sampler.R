sbm.gibbs <- function(A,z.init,th.init,mu.init,n.iters,K,mu0,sig0,alpha){
    # setup a Gibbs sampler for weighted SBM
    # A -       n x n adjacency matrix
    # z.init -  initial values for z (vertex labels)
    # th.init - K x K initial values for theta (block parameters)
    # mu.init - initial values for mu (community assignment probabilities)
    # n.iters - how many draws to take
    # K -       number of communities
    # mu0 -     mu_0 parameter of Normal prior for theta
    # sig0 -    sig_0^2 parameter of Normal prior for theta
    # alpha -   alpha parameter of Dirichlet prior for mu
    # for now, hardcode to use normal prior, Dirichlet prior


    n <- dim(A)[1]
    
    ## initialize to hold iterations
    z.p <- array(0,dim=c(n.iters,n,K)) # iters by n by K (1-hots)
    
    th.p <- array(0,dim=c(n.iters,K,K)) # iters by K by K

    mu.p <- matrix(0,nrow=n.iters,ncol=length(mu.init)) # iters by K


    draw.params <- function(p){
        ## draw the parameters conditionally on everything else
        params <- list()
        params$th <- draw.th(z.p[p,,])
        params$mu <- draw.mu(z.p[p,,])

        return(params)

    }

    draw.th <- function(z){
        ## draw theta conditionally on everything else
        ns <- colSums(z) # how many nodes in each community
        bns <- ns %o% ns # how many edges in each block
        ts <- t(z) %*% A %*% z # sum edges in each block; seems risky

        mus <- 1/(1/sig0 + bns/1) * (mu0/sig0 + ts/1) # posterior mus
        sigs <- 1/(1/sig0 + bns/1) # posterior sigma^2

        th <- matrix(rnorm(K^2,mus,sigs),nrow=K,ncol=K)

        return(th)
    }

    rdir <- function(alpha){
        ## sample from the Dirichlet distribution, not in base R
        n <- length(alpha)
        gams <- rgamma(n,alpha)
        x <- gams/sum(gams)
        return(x)
    }

    draw.mu <- function(z){
        ## draw mu conditionally on everything else
        ns <- colSums(z)
        alphas <- alpha + ns # posterior alpha
        mu <- rdir(alphas)
        return(mu)
    }



    A.ll <- function(z,th){
        ## return the log likelihood of A given z, th
        A.ctr <- A - z %*% th %*% t(z) # center A about block-specific 
        ll <- t(z) %*% (-A.ctr^2/2) %*% z    # sum over negate square centered A divided 2 sig^2
        l <- sum(ll)                         # sum over the logs
        return(l)
        }
        

    draw.z <- function(p){
        ## draw all the z's conditionally with a nested gibbs thing

        z.new <- matrix(nrow=n,ncol=K) # we'll send this back

        z.cur <- z.p[p-1,,] # first one is special
        z.cur[1,] <- 0
        z.new[1,] <- draw.z1(1,z.cur)

        for (i in 2:(n-1)){                                # intermediate ones
            z.cur <- rbind(z.new[1:i,], z.p[p-1,(i+1):n,]) # use any new Z's we can
            z.cur[i,] <- 0                                 # make sure it's 0
            z.new[i,] <- draw.z1(i,z.cur)                  # get the new value
        }

        z.cur <- z.new[,] # last one is special
        z.cur[n,] <- 0
        z.new[n,] <- draw.z1(n,z.cur)

        return(z.new)
    }

    draw.z1 <- function(i,z.cur){
        ## draw the i'th z conditionally
        probs <- rep(0,K) # initialize
        for (k in 1:K){   # proportional to each possible value
            z.try <- z.cur
            z.try[i,k] <- 1
            probs[k] <- A.ll(z.try,th.p[p,,])
        }
        probs <- probs/sum(probs) # push it onto simplex
        zi <- rmultinom(1,1,probs)
        return(zi)
    }

    constrain.space <- function(th,mu,z){
        ## make sure that the diagonal of th is descending
        idx <- order(diag(th))

        new.th <- th[idx,idx]
        new.mu <- mu[idx]
        new.z <- z[,idx]

        return(list(th=new.th,mu=new.mu,z=new.z))
        }
    
    ## seed the first iterate
    z.p[1,,] <- z.init
    th.p[1,,] <- th.init
    mu.p[1,] <- mu.init

    p <- 2

    for (p in 2:n.iters){

        params.new <- draw.params(p)
        th.p[p,,] <- params.new$th
        mu.p[p,] <- params.new$mu

        z.new <- draw.z(p)
        z.p[p,,] <- z.new

        ## permute stuff to make diagonal descending
        ordered <- constrain.space(th.p[p,,],mu.p[p,],z.p[p,,])
        th.p[p,,] <- ordered$th
        mu.p[p,] <- ordered$mu
        z.p[p,,] <- ordered$z
    }



    return(list(z.p=z.p,th.p=th.p,mu.p=mu.p))
}
