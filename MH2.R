# Save initial time
time0 <- Sys.time()

load("13253_2017_282_MOESM2_ESM.rdata")
source("MH2_fun.R")
source("mllk.R")

# Retrieve initial parameter values from "mod"
load("par.vec0.RData")
pn <- mod$pn
npar <- length(as.vector(unlist(mod$pn)))
pn <- rapply(pn, function(x) round(x, digits=1), how="replace") # approximate pn

n.iter <- 16000 # iterations
n.tune <- 6000 # iterations for tuning

# Set sd of the normal proposal distribution for random walk MH update
delta0 <- matrix(c(NA, NA, NA, NA, 0.18, 0.16, 0.10, 0.07, 0.05, 0.06, 0.02, NA, NA, NA, NA, 0.60,
                   0.44, 0.26, 0.19, 0.30, 0.32, 0.003, NA, NA, NA, NA, 2.00, 1.00, 0.73, 0.40, 
                   1.10, 0.67, 0.001), 
                 ncol = 3, byrow=F)

# Set thresholds for acceptance rates
thr <- c(0.25, 0.4)

# Parameter values for the priors of dw.pi
alpha <- rep(1, 3)
beta <- rep(1, 3)

# Parameter values for the priors of means
mu  <- rep(log(100), 3)
sd <- rep(1, 3)

# Parameter values for the priors of sd's
shape  <- rep(1e-3, 3)
rate <- rep(1e-3, 3)

# Calculate log-likelihood for initial state 
likelihood <- -mllk(pn=pn, data=data, ll.N=2, ul.N=3, fit=TRUE)

# Perform MH algorithm with block updates
output <- updatepar_MH2(par.list=pn, npar=npar, data=data, likhood=likelihood, 
                      alpha=alpha, beta=beta, delta=delta0, mu=mu, sigma=sd,
                      shape=shape, rate=rate, n.iter, tuning=T, thr=thr, n.tune=n.tune)
  
# Save output
itns_MH2 <- output$itns[31:51,]
acc.count <- output$counts
delta <- output$delta

# Proportion of accepted moves over iterations
acc.rate <- apply(acc.count, c(1,2), sum)/n.iter

# Monitor computational time
timef <- Sys.time()

# Save estimates
estimates <- rowMeans(itns_MH2[,(n.iter-n.tune+1):n.iter])
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
