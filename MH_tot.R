load("13253_2017_282_MOESM2_ESM.rdata")
source("updatepar.R")
source("mllk.R") 

# Load initial parameter values from "mod"
load("par.vec0.RData") 
pn <- mod$pn
npar <- length(as.vector(unlist(mod$pn)))
pn <- rapply(pn, function(x) round(x, digits=3), how="replace") ## Approximate pn

n.iter <- 16000 ## Total iterations
n.tune <- 6000 ## Pilot tuning

# Set sd of the normal proposal distribution for random walk MH update
delta <- list(ll.delta=c(0.05,0.05),
              ll.gamma=matrix(rep(0.03,4), nrow=2, byrow=T),
              ul.delta=list(rep(0.02,3), rep(0.02,3)),
              ul.gamma=list(matrix(c(0.02,0.02,0.02,
                                     0.02,0.02,0.02, 
                                     0.01,0.01,0.01), byrow=T, nrow=3),
                            matrix(c(0.013,0.013,0.013, 
                                     0.013,0.013,0.013,
                                     0.01,0.01,0.01), byrow=T, nrow=3)),
              dd.mu=c(0.2,0.60,2.00),
              dd.sigma=c(0.15,0.5,1.0),
              md.mu=c(0.10,0.25,0.70),
              md.sigma=c(0.08,0.20,0.5),
              dw.mu=c(0.05,0.30,1.10),
              dw.sigma=c(0.06,0.30,0.7),
              dw.pi=c(0.02,0.003, 0.001))

# Create a list to store the count of accepted move for each parameter through all the iterations
# the number of iterations will be updated at each iteration as a sum
acc.count <- rapply(pn, function(x)ifelse(x!=0,0,x), how="replace")

# Parameter values for the priors of dw.pi
alpha <- rep(1, 3)
beta <- rep(1, 3)

# Parameter values for the priors of means
mu  <- rep(log(100), 3)
sd <- rep(1, 3)

# Parameter values for the priors of sd's
shape  <- rep(1e-3, 3)
rate <- rep(1e-3, 3)

thr <- c(0.25, 0.40)

time0 <- Sys.time() ## Save initial time

# Calculate log-likelihood for initial state 
likelihood <- -mllk(pn=pn, data=data, ll.N=2, ul.N=3, fit=TRUE)

# Define an array to store sample from posterior distribution for each chain
itns <- array(0, dim=c(npar,n.iter))
row.names(itns) <- names(unlist(pn))

# MCMC updates - MH algorithm - cycle through each iteration
for (t in 1:n.iter){
  
  output <- updatepar_tot(par.list=pn, npar=npar, data=data, likhood=likelihood, alpha=alpha, 
                          beta=beta, delta=delta, mu=mu, sigma=sd, shape=shape, rate=rate,
                          acc.count=acc.count)
  
  pn <- output$par.list 
  likelihood <- output$likhood
  acc.count <- output$acc.count
  
  itns[,t] <- as.vector(unlist(pn))
  
  if (t<n.tune){ ## tuning
    if (mod(t,100)==0){ ## every 100 iterations

      foo <- as.vector(unlist(acc.count))/t
      delta <- as.vector(unlist(delta))

      delta_tmp <- ifelse(foo < thr[1], 0.9*delta,
                          ifelse(foo > thr[2], 1.1*delta, delta))

      delta <- list(ll.delta=delta[1:2],
                    ll.gamma=matrix(delta[3:6], nrow=2),
                    ul.delta=list(delta[7:9],delta[10:12]),
                    ul.gamma=list(matrix(delta[13:21],nrow=3),
                                  matrix(delta[22:30],nrow=3)),
                    dd.mu=delta_tmp[31:33],
                    dd.sigma=delta_tmp[34:36],
                    md.mu=delta_tmp[37:39],
                    md.sigma=delta_tmp[40:42],
                    dw.mu=delta_tmp[43:45],
                    dw.sigma=delta_tmp[46:48],
                    dw.pi=delta_tmp[49:51])
    }
  }
}

timef <- Sys.time() ## Monitor computational time
tot_time <- time0 - timef

estimates <- rowMeans(itns[,(n.iter-n.tune+1):n.iter])
est <- list(ll.delta=estimates[1:2],
            ll.gamma=matrix(estimates[3:6],ncol=2),
            ul.delta=list(estimates[7:9], estimates[10:12]),
            ul.gamma=list(matrix(estimates[13:21],nrow=3),
                          matrix(estimates[22:30],nrow=3)),
            dd.mu=estimates[31:33],
            dd.sigma=estimates[34:36],
            md.mu=estimates[37:39],
            md.sigma=estimates[40:42],
            dw.mu=estimates[43:45],
            dw.sigma=estimates[46:48],
            dw.pi=estimates[49:51])
