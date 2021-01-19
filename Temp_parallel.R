library(foreach)
library(doParallel)

time0 <- Sys.time() # Save initial time

# Build clusters
no_cores <- detectCores() - 1
#no_cores <- 4 # use 4 cores
cl <- makeCluster(no_cores)
registerDoParallel(cl)

load("13253_2017_282_MOESM2_ESM.rdata") # Load data
source("Temp_parallel_fun.R") # Load likelihood function

# Initialize emperature values for the chains (for all parameters)
B.temp <- seq(1, 0, length.out=(no_cores+1))[(-(no_cores+1))]
#B.temp <- c(1, 0.75, 0.5, 0.25) # use this for 4 cores
n.chains <- length(B.temp)

# Retrieve initial parameter values from "mod"
load("par.vec0.RData")
pn <- mod$pn
n.par <- length(as.vector(unlist(mod$pn)))
pn <- rapply(pn, function(x) round(x, digits=1), how="replace") # Approximate pn

n.iter <- 100 # set the number of iterations (between swaps)
n.swaps <- 160 # set the number of swaps
n.tune <- 6000 # set the length of the tuning

# Parameter values for the priors of dw.pi
alpha <- rep(1, 3)
beta <- rep(1, 3)

# Parameter values for the priors of means
mu  <- rep(log(100), 3)
sd <- rep(1, 3)

# Parameter values for the priors of sd's
shape  <- rep(1e-3, 3)
rate <- rep(1e-3, 3)

# Define a list of sd's for the RW 
delta <- matrix(c(NA, NA, NA, NA, 0.18, 0.16, 0.10, 0.07, 0.05, 0.06, 0.02, NA, NA, NA, NA, 0.60,
                  0.44, 0.26, 0.19, 0.30, 0.32, 0.003, NA, NA, NA, NA, 2.00, 1.00, 0.73, 0.40, 
                  1.10, 0.67, 0.001), ncol = 3, byrow=F)
colnames(delta) <- cols <-cols <- c("ps1", "ps2","ps3")
rownames(delta) <- rows <- c("ll.delta", "ll.gamma", "ul.delta", "ul.gamma", "dd.mu",
                             "dd.sigma", "md.mu", "md.sigma", "dw.mu", "dw.sigma", "dw.pi")

# Set thresholds for acceptance rates
thr <- c(0.25, 0.4)

# Calculate log-likelihood for initial state 
likelihood <- rep(-mllk(pn=mod$pn, data=data, ll.N=2, ul.N=3, fit=TRUE),n.chains)

# Initialize the value of the energy (and keep it as a temporary storage)
E <- -likelihood - prior(par.list=mod$pn, mu=mu, sigma=sd, alpha=alpha, beta=alpha, shape=shape, rate=rate)

# Create a dataset to store the count of accepted move for each parameter through all the iterations
# the number of iterations will be updated at each iteration as a sum
# rows=parameters, cols=prod states
acc.count <- array(0, dim=c(11,3, n.iter))

# Define an array to store sample from posterior distribution for the chosen "final" chain
itns_pt <- matrix(0, nrow=n.par, ncol=(n.iter*n.swaps))
# Array of iterations for each chain
itns.tmp <- array(0, dim=c(n.par, n.iter))

# Repeat for all the swaps
for (s in 1:n.swaps){
  
  result <- foreach(j=1:n.chains,
                    .packages=c("RcppArmadillo", "Rcpp", "boot", "invgamma", "numbers"),   # instead of clusterEvalQ 
                    .combine=list,
                    .multicombine=TRUE) %do%
    temp.tot.update(j)
  
  for (j in n.chains:2){
    
    # Pairwise difference of temperature values for the chains
    D.B <- abs(result[[j]]$B.temp - result[[j-1]]$B.temp)
    
    # Pairwise differences in energy
    D.E <- as.vector(result[[j]]$energy - result[[j-1]]$energy)
    
    # Acceptance prob for the swap for each parameter
    A <- min(1, exp(-D.B*D.E))
    
    if(runif(1)<= A){  # accept the chains swap
      
      # swap B, E, itns, likelihood, pn dvalues accordingly to the chain swap
      B.tmp <- result[[j]]$B.temp
      result[[j]]$B.temp <- result[[j-1]]$B.temp
      result[[j-1]]$B.temp <- B.tmp
      
      E.tmp <- result[[j]]$energy
      result[[j]]$energy <- result[[j-1]]$energy
      result[[j-1]]$energy <- E.tmp
      
      L.tmp <- result[[j]]$likelihood
      result[[j]]$likelihood <- result[[j-1]]$likelihood
      result[[j-1]]$likelihood <- L.tmp
      
      I.tmp <- result[[j]]$itns
      result[[j]]$itns <- result[[j-1]]$itns
      result[[j-1]]$itns <- I.tmp
      
      pn.tmp <- result[[j]]$pn
      result[[j]]$pn <- result[[j-1]]$pn
      result[[j-1]]$pn <- pn.tmp
      
    }
  }
  
  # Save pn and likelihood for the selected chain
  pn <- result[[1]]$pn
  likelihood <- result[[1]]$likelihood
  
  # Update itns based on the values of the selected chain
  itns_pt[,((s-1)*n.iter+1) : (s*n.iter)] <- result[[1]]$itns
  
  # Save the vector of temperatures and energies for all the chains
  B.temp <- sapply(1:n.chains, function(x){result[[x]]$B.temp})
  E <- sapply(1:n.chains, function(x){result[[x]]$energy})
  
  clusterExport(cl, list("pn", "B.temp", "E", "likelihood"))
  
}

timef <- Sys.time() # Monitor computational time
stopCluster(cl) # Close clusters

# Save estimates
estimates <- rowMeans(itns_pt[,(n.iter*n.swaps-n.tune+1):(n.iter*n.swaps)])
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

