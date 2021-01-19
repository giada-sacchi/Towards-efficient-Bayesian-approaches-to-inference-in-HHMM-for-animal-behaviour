library(foreach)
library(doParallel)

time0 <- Sys.time()

# Build clusters
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

load("13253_2017_282_MOESM2_ESM.rdata")
source("Temp_parallel_fun_TOT.R")

# Initialize emperature values for the chains (for all parameters)
B.temp <- c(1, 0.75, 0.5, 0.25)
n.chains <- length(B.temp)

# Retrieve initial parameter values from "mod"
load("par.vec0.RData")
pn <- mod$pn
n.par <- length(as.vector(unlist(mod$pn)))
pn <- rapply(pn, function(x) round(x, digits=3), how="replace")

n.iter <- 100 # Set the number of iterations (between swaps)
n.swaps <- 16 # set the number of swaps
n.tune <- 600 # set the length of the tuning

# Parameter values for the priors of dw.pi
alpha <- rep(1, 3)
beta <- rep(1, 3)

# Parameter values for the priors of means
mu  <- rep(log(100), 3)
sd <- rep(1, 3)

# Parameter values for the priors of sd's
shape  <- rep(1e-3, 3)
rate <- rep(1e-3, 3)

# define a list of sd's for the RW 
delta <- list(ll.delta=c(0.05,0.05),
              ll.gamma=matrix(rep(0.03,4), nrow=2, byrow=T),
              ul.delta=list(rep(0.02,3), rep(0.02,3)),
              ul.gamma=list(matrix(c(0.02,0.02,0.02,
                                     0.02,0.02,0.02, 
                                     0.01,0.01,0.01),byrow=T,nrow=3),
                            matrix(c(0.013,0.013,0.013, 
                                     0.013,0.013,0.013,
                                     0.01,0.01,0.01),byrow=T,nrow=3)),
              dd.mu=c(0.2,0.60,2.00),
              dd.sigma=c(0.15,0.5,1.0),
              md.mu=c(0.10,0.25,0.70),
              md.sigma=c(0.08,0.20,0.5),
              dw.mu=c(0.05,0.30,1.10),
              dw.sigma=c(0.06,0.30,0.7),
              dw.pi=c(0.02,0.003, 0.001))

# Set thresholds for acceptance rates
thr <- c(0.25, 0.4)

# Calculate log-likelihood for initial state 
likelihood <- -mllk(pn=pn, data=data, ll.N=2, ul.N=3, fit=TRUE)

# Initialize the value of the energy (and keep it as a temporary storage)
E <- -likelihood - prior(par.list=pn, mu=mu, sigma=sd, alpha=alpha, beta=alpha, shape=shape, rate=rate)

# Create a list to store the count of accepted move for each parameter through all the iterations
# the number of iterations will be updated at each iteration as a sum
acc.count <- rapply(pn, function(x)ifelse(x!=0,0,x), how="replace")

# Define an array to store sample from posterior distribution for the chosen "final" chain
itns_pt <- matrix(0, nrow=n.par, ncol=(n.iter*n.swaps))
row.names(itns_pt) <- names(unlist(pn))
# Array of iterations for each chain
itns.tmp <- array(0, dim=c(n.par, n.iter))

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
      
      # swap B, E, itns, likelihood values accordingly to the chain swap
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
      
      counts.tmp <- result[[j]]$counts
      result[[j]]$counts <- result[[j-1]]$counts
      result[[j-1]]$counts <- counts.tmp
    }
  }
  
  pn <- result[[1]]$pn
  B.temp <- sapply(1:n.chains, function(x){result[[x]]$B.temp})
  E <- sapply(1:n.chains, function(x){result[[x]]$energy})
  likelihood <- result[[1]]$likelihood
  itns_pt[,((s-1)*n.iter+1):(s*n.iter)] <- result[[1]]$itns
  acc.count <- result[[1]]$counts
  
  clusterExport(cl, list("pn", "B.temp", "E", "likelihood","acc.count"))
  
}

timef <- Sys.time()
stopCluster(cl)

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

