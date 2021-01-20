time0 <- Sys.time()

load("13253_2017_282_MOESM2_ESM.rdata")
source("MH-algo_temp_fun.R")
source("mllk.R")

# Retrieve initial parameter values from "mod"
load("par.vec0.RData")
pn <- mod$pn
n.par <- length(as.vector(unlist(mod$pn)))
pn <- rapply(pn, function(x) round(x, digits=1), how="replace") # approximate pn
pn <- list(pn, pn, pn, pn)

n.iter <- 100 # iterations between swaps
n.swaps <- 100 # number of swaps -> 10'000 total iterations
n.tune <- 6000 # iterations out of n.iter*n.swaps on which performing pilot tuning

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
delta <- matrix(c(NA, NA, NA, NA, 0.18, 0.16, 0.10, 0.07, 0.05, 0.06, 0.02, NA, NA, NA, NA, 0.60,
                  0.44, 0.26, 0.19, 0.30, 0.32, 0.003, NA, NA, NA, NA, 2.00, 1.00, 0.73, 0.40, 
                  1.10, 0.67, 0.001), ncol = 3, byrow=F)
colnames(delta) <- cols <- c("ps1", "ps2","ps3")
rownames(delta) <- rows <- c("ll.delta", "ll.gamma", "ul.delta", "ul.gamma", "dd.mu",
                             "dd.sigma", "md.mu", "md.sigma", "dw.mu", "dw.sigma", "dw.pi")

# Set thresholds for acceptance rates
thr <- c(0.25, 0.4)

# Initialize emperature values for the chains (for all parameters)
B.temp <- c(1, 0.75, 0.5, 0.25)
n.chains <- length(B.temp)

# Calculate log-likelihood for initial state (one for each chain)
likelihood <- rep(-mllk(pn=pn[[1]], data=data, ll.N=2, ul.N=3, fit=TRUE), n.chains)

# Create a dataset to store the count of accepted move for each parameter through all the iterations
# the number of iterations will be updated at each iteration as a sum
# rows=parameters, cols=prod states
acc.count <- array(0, dim=c(11,3, n.iter))

# Define an array to store sample from posterior distribution for each chain
itns <- matrix(0, nrow=n.par, ncol=(n.iter*n.swaps))
itns.tmp <- array(0, dim=c(n.par, n.iter, n.chains))

# Initialize the value of the energy (and keep it as a temporary storage)
E <- -likelihood - prior(par.list=pn[[1]], mu=mu, sigma=sd, alpha=alpha, beta=beta, shape=shape, rate=rate)

# MCMC updates - MH algorithm
for (s in 1:n.swaps){
  
  for (j in 1:n.chains){
    
    for (t in 1:n.iter){

      output <- updatepar.temp(par.list=pn[[1]], npar=n.par, data=data, likhood=likelihood[1], 
                               shape=shape, rate=rate, alpha=alpha, beta=beta, delta=delta,
                               mu=mu, sigma=sd, B.temp=B.temp[j], E=E[j])
      
      pn[[j]] <- output$par.list 
      likelihood[j] <- output$likhood
      acc.count[,,t] <- output$counts
      E[j] <- output$energy      # save the final value of energy (overwrites until the (s+1)-th swap)

      # save temporarily iterations for each chain
      itns.tmp[,t,j] <- as.vector(unlist(pn[[j]]))
      
      if(t<n.tune){
        if (mod(t,100)==0){ # every 100 iterations
          
          #look at accepted moves for the previous t iterations
          foo <- apply(acc.count, c(1,2), sum)/t
          
          # if less than thr[1]% of the moves were accepted, reduce delta (only for that parameter)
          delta[foo < thr[1]] <- 0.9*delta[foo < thr[1]]
          
          # if more than thr[2]% of the moves were accepted, increase delta (only for that parameter)
          delta[foo > thr[2]] <- 1.1*delta[foo > thr[2]]
        }
      }
    }
  }
  
  for (j in n.chains:2){
    
    # Pairwise difference of temperature values for the chains
    D.B <- abs(B.temp[j] - B.temp[j-1])
    
    # Pairwise differences in energy
    D.E <- as.vector(E[j] - E[j-1])
    
    # Acceptance prob for the swap for each parameter
    A <- min(1, exp(-D.B*D.E))
    
    if(runif(1)<= A){  # accept the chains swap
      
      # swap B, E, itns, likelihood values accordingly to the chain swap
      B.temp[c(j-1,j)] <- B.temp[c(j,j-1)]
      E[c(j-1,j)] <- E[c(j,j-1)]
      likelihood[c(j-1,j)] <- likelihood[c(j,j-1)]
      itns.tmp[,,c(j,j-1)] <- itns.tmp[,,c(j-1,j)]
      
      pn.tmp <- pn[[j]]
      pn[[j]] <- pn[[j-1]]
      pn[[j-1]] <- pn.tmp
      
    }
  }
  
  # just save new iterations (at "t") for the chosen chain
  itns[,((s-1)*n.iter+1) : (s*n.iter)] <- itns.tmp[, , 1]
  
}

itns_temp <- itns
rm(itns)
timef <- Sys.time()

estimates <- rowMeans(itns_temp[,(n.iter*n.swaps-n.tune):(n.iter*n.swaps)])
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
