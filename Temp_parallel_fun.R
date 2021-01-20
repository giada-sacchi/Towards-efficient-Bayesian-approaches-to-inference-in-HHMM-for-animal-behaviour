##### UPDATEPAR.TEMP #####
# Function used to update parameters (acceptance prob computed considering chain temperature)
updatepar.temp <- function(par.list, npar, data, 
                      likhood,                     # current likelihood
                      B.temp=1,                    # temperature of the chain
                      E,                           # initial energy of the chain
                      mu, sigma, alpha, beta, shape, rate,
                      delta){            # matrix with sd of random walk proposal distribution 
  # counts of the accepted moves for each parameter (rows=parameters, cols=prod states)
  counts <- matrix(rep(0,33), ncol=3)
  for (j in 5:11){ # positions corresponding to target variable parameters
    for (i in 1:3) { # for the three prod states
      # Keep a record of the current parameter value being updated 
      oldpar <- par.list[[j]][i]
      # Propose a new candidate value using a random walk from normal proposal 
      par.list[[j]][i] <- rnorm(1, par.list[[j]][i], delta[j,i])
      if (j==11) { # dw.pi
        if (par.list[[j]][i]>=0 & par.list[[j]][i]<=1){ # dw.pi is a probability
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
            # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
            num <- newlikhood + dbeta(par.list[[j]][i], alpha[i], beta[i],log=TRUE)
            den <- likhood + dbeta(oldpar, alpha[i], beta[i],log=TRUE)
            E <- den-num
            A <- min(1, exp(-B.temp*E))
          }
        }  else { # if parameter is outside range set the acceptance probability A to 0
          A <- 0 
        }
        A
      } else if (is.element(j,c(5,7,9))){ # parameter is a mean
        if (par.list[[j]][i]>=0){
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
          # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
          num <- newlikhood + dlnorm(par.list[[j]][i], mu[i], sd[i],log=TRUE)
          den <- likhood + dlnorm(oldpar, mu[i], sd[i],log=TRUE)
          E <- den-num
          A <- min(1, exp(-B.temp*E))
          }
        }  else { # Otherwise set the acceptance probability A to 0
          A <- 0 
        }
      } else { # parameter is a sd
        if (par.list[[j]][i]>=0){
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
          # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
          num <- newlikhood + dinvgamma(par.list[[j]][i], shape[i], rate[i],log=TRUE)
          den <- likhood + dinvgamma(oldpar, shape[i], rate[i],log=TRUE)
          E <- den-num
          A <- min(1, exp(-B.temp*E))
          }
        }  else { # Otherwise set the acceptance probability A to 0
          A <- 0 
        }
      }
      if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
        # Accept the move with probability A and store its likelihood
        likhood <- newlikhood
        counts[j,i] <- 1
      } else {     
        # Reject move and return parameter value to previous value
        par.list[[j]][i] <- oldpar
      }
    }
  }
  # Output the parameter values and log-likelihood value (and the vector of counts for accepted moves)
  output <- list(par.list=par.list, likhood=likhood, energy=E, counts=counts) 
  output
}


##### PRIOR #####
# Function to compute the sum of log-prior distributions of a parameter vector

prior <- function(par.list, mu, sigma, alpha, beta, shape, rate){
  prior <- 0
  for (j in 5:11){ # positions corresponding to target variable parameters
    for (i in 1:3) {
      if (j==11) { 
        prior <- prior + dbeta(par.list[[j]][i], alpha[i], beta[i],log=TRUE)
        } else if (is.element(j, c(5,7,9))){
          prior <- prior + dlnorm(par.list[[j]][i], mu[i], sigma[i],log=TRUE)
        } else {  
          prior <- prior + dinvgamma(par.list[[j]][i], shape[i], rate[i],log=TRUE)
        }
    }
  }
  return(prior)
}

##### TEMP.TOT.UPDATE #####
# create function for updating the parameter vector "parvec" over "n" iterations
temp.tot.update <- function(j){
  sourceCpp("13253_2017_282_MOESM4_ESM.cpp")
  B <- B.temp[j]
  E <- E[j]
  # parameter update
  for (t in 1:n.iter){
    output <- updatepar.temp(par.list=pn, npar=n.par, data=data, likhood=likelihood, 
                             shape=shape, rate=rate, alpha=alpha, beta=beta, delta=delta,
                             mu=mu, sigma=sd, B.temp=B, E=E)
    pn <- output$par.list 
    likelihood <- output$likhood
    acc.count[,,t] <- output$counts
    E <- output$energy
    itns.tmp[,t] <- as.vector(unlist(pn))  # save temporarily iterations for each chain
    if(t<n.tune){ # tuning
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
  return(list(pn=pn, likelihood=likelihood, energy=E, itns=itns.tmp, B.temp=B))
}
