load("13253_2017_282_MOESM2_ESM.rdata")
source("MH_fun.R")
source("mllk.R")

# Save initial time
time0 <- Sys.time()

# Retrieve initial parameter values from "mod"
load("par.vec0.RData")
pn <- mod$pn
npar <- length(as.vector(unlist(mod$pn)))
pn <- rapply(pn, function(x) round(x, digits=1), how="replace") # Approximate pn

# Set the number of iterations
n.iter <- 16000

# Set number of iterations for tuning
n.tune <- 6000

# For the first trial
if (!exists("acc.rate")){
  
  # Initiate an empty list in which store datasets for counting accepted moves for each parameters
  # Each sub-list will correspond to a different set of delta
  acc.rate <- vector("list", 1)
  
  # Set sd of the normal proposal distribution for random walk MH update
  delta <- matrix(c(NA, NA, NA, NA, 0.18, 0.16, 0.10, 0.07, 0.05, 0.06, 0.02, NA, NA, NA, NA, 0.60,
                    0.44, 0.26, 0.19, 0.30, 0.32, 0.003, NA, NA, NA, NA, 2.00, 1.00, 0.73, 0.40, 
                    1.10, 0.67, 0.001), 
                  ncol = 3, byrow=F)
  
  acc.moves <- matrix(rep(0, 33), nrow=11)
  colnames(acc.moves) <- c("ps1","ps2","ps3")
  rownames(acc.moves) <- c("ll.delta", "ll.gamma", "ul.delta", "ul.gamma", "dd.mu",
                           "dd.sigma", "md.mu", "md.sigma", "dw.mu", "dw.sigma", "dw.pi")
  
  # Store the same value in the list for recording counts and a vector with the count of accepted moves
  acc.rate[[length(acc.rate)]] <- list(delta=delta, acc.moves=acc.moves, thr=c(0.25, 0.4))  
}

# Create a dataset to store the count of accepted move for each parameter through all the iterations
# the number of iterations will be updated at each iteration as a sum
# rows=parameters, cols=production states
acc.count <- array(0, dim=c(11,3, n.iter))

# Parameter values for the priors of dw.pi
alpha <- rep(1, 3)
beta <- rep(1, 3)

# Parameter values for the priors of means
mu  <- rep(log(100), 3)
sd <- rep(1, 3)

# Parameter values for the priors of sd's
shape  <- rep(1e-3, 3)
rate <- rep(1e-3, 3)

# Set thresholds for delta if desired different from default c(0.3,0.5)
#acc.rate[[length(acc.rate)]]$thr <- thr <- c(0.25, 0.40)

# Calculate log-likelihood for initial state 
likelihood <- -mllk(pn=pn, data=data, ll.N=2, ul.N=3, fit=TRUE)

# Define an array to store sample from posterior distribution for each chain
itns <- array(0, dim=c(n.iter, npar))

# MCMC updates - MH algorithm - cycle through each iteration
for (t in 1:n.iter){
  
  output <- updatepar_MH(par.list=pn, npar=npar, data=data, likhood=likelihood, alpha=alpha, beta=beta, 
                         delta=delta, mu=mu, sigma=sd, shape=shape, rate=rate)
  # Save output
  pn <- output$par.list 
  likelihood <- output$likhood
  acc.count[,,t] <- output$counts
  
  itns[t,] <- as.vector(unlist(pn))
  
  # perform tuning
  if (t<n.tune){
    if (mod(t,100)==0){ # every 100 iterations
      
      thr <- acc.rate[[length(acc.rate)]]$thr
      
      #look at accepted moves for the previous t iterations
      foo <- apply(acc.count, c(1,2), sum)/t
      
      # if less than thr[1]% of the moves were accepted, reduce delta (only for that parameter)
      delta[foo < thr[1]] <- 0.9*delta[foo < thr[1]]
      
      # if more than thr[2]% of the moves were accepted, increase delta (only for that parameter)
      delta[foo > thr[2]] <- 1.1*delta[foo > thr[2]]
    }
  }
}

# Update the count of accepted moves per parameter for the current value of delta
acc.rate[[length(acc.rate)]]$acc.moves <- apply(acc.count, c(1,2), sum)/n.iter

# Create a new list for the next value of delta
par.names <- c("ll.delta", "ll.gamma", "ul.delta", "ul.gamma", "dd.mu", "dd.sigma", "md.mu",
               "md.sigma", "dw.mu", "dw.sigma", "dw.pi")
acc.rate[[length(acc.rate)+1]] <- list(delta=delta, 
                                       acc.moves=matrix(rep(0, 33), nrow=11, 
                                                        dimnames=list(par.names, c("ps1","ps2","ps3"))))
# Monitor computational time
timef <- Sys.time()

# Save estimates
estimates <- colMeans(itns[(n.iter-n.tune+1):n.iter,])
est <- list(ll.delta=estimates[1:2],
            ll.gamma=estimates[3:6],
            ul.delta=list(estimates[7:9], estimates[10:12]),
            ul.gamma=list(estimates[13:21],estimates[22:30]),
            dd.mu=estimates[31:33],
            dd.sigma=estimates[34:36],
            md.mu=estimates[37:39],
            md.sigma=estimates[40:42],
            dw.mu=estimates[43:45],
            dw.sigma=estimates[46:48],
            dw.pi=estimates[49:51])
