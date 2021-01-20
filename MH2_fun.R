##### UPDATEPAR_MH2 #####
# Function for single-update algorithm (based on Block2)
updatepar_MH2 <- function(par.list, npar, data, n.iter, tuning=F, thr=c(0.2,0.4),
                             likhood,       # current likelihood
                             mu, sigma, alpha, beta, shape, rate, n.tune=n.iter,
                             delta){        # matrix with sd of random walk proposal distribution 
  # Define an array to store sample from posterior distribution for each chain
  itns <- array(0, dim=c(npar, n.iter))
  # counts of the accepted moves for each parameter (rows=parameters, cols=prod states)
  counts <- array(0, dim=c(11,3, n.iter))
  for (t in 1:n.iter){
    for (i in 1:3) { # for the three prod states
      # Keep a record of the current parameter value being updated 
      oldpar <- par.list[[11]][i]
      # Propose a new candidate value using a random walk from normal proposal 
      par.list[[11]][i] <- rnorm(1, par.list[[11]][i], delta[11,i])
      if (par.list[[11]][i]>=0 & par.list[[11]][i]<=1){ # dw.pi is a probability
        # Compute the log-likelihood (of the move)
        newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
        if(is.nan(newlikhood)){
          A <- 0
        } else {
          # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
          num <- newlikhood + log(dbeta(par.list[[11]][i], alpha[i], beta[i]))
          den <- likhood + log(dbeta(oldpar, alpha[i], beta[i]))
          A <- min(1, exp(num-den))
        }
      }  else { # if parameter is outside range set the acceptance probability A to 0
        A <- 0 
      }
      if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
        # Accept the move with probability A and store its likelihood
        likhood <- newlikhood
        counts[11,i,t] <- 1
      } else {     
        # Reject move and return parameter value to previous value
        par.list[[11]][i] <- oldpar
      }
    }
    for (j in c(5,7,9)){ # parameter is a mean
      for (i in 1:3) { # for the three prod states
        # Keep a record of the current parameter value being updated 
        oldpar <- par.list[[j]][i]
        # Propose a new candidate value using a random walk from normal proposal 
        par.list[[j]][i] <- rnorm(1, par.list[[j]][i], delta[j,i])
        if (par.list[[j]][i]>=0){
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
          num <- newlikhood + log(dlnorm(par.list[[j]][i], mu[i], sd[i]))
          den <- likhood + log(dlnorm(oldpar, mu[i], sd[i]))
          A <- min(1, exp(num-den))
        }  else { # Otherwise set the acceptance probability A to 0
          A <- 0 
        }
        if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
          # Accept the move with probability A and store its likelihood
          likhood <- newlikhood
          counts[j,i,t] <- 1
        } else {     
          # Reject move and return parameter value to previous value
          par.list[[j]][i] <- oldpar
        }
      }
    }
   for (j in c(6,8,10)) { # parameter is a sd
      for (i in 1:3) { # for the three prod states
        # Keep a record of the current parameter value being updated 
        oldpar <- par.list[[j]][i]
        # Propose a new candidate value using a random walk from normal proposal 
        par.list[[j]][i] <- rnorm(1, par.list[[j]][i], delta[j,i])
        if (par.list[[j]][i]>=0){
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
            # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
            num <- newlikhood + log(dinvgamma(par.list[[j]][i], shape[i], rate[i]))
            den <- likhood + log(dinvgamma(oldpar, shape[i], rate[i]))
            A <- min(1, exp(num-den))
          }
        }  else { # Otherwise set the acceptance probability A to 0
          A <- 0 
        }
        if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
          # Accept the move with probability A and store its likelihood
          likhood <- newlikhood
          counts[j,i,t] <- 1
        } else {     
          # Reject move and return parameter value to previous value
          par.list[[j]][i] <- oldpar
        }
      }
    }
    itns[,t] <- as.vector(unlist(par.list))
    if (tuning==T){
      if (t<n.tune){
        if (mod(t,100)==0){ # every 100 iterations
          #look at accepted moves for the previous t iterations
          foo <- (apply(counts, c(1,2), sum)/t)
          # if less than thr[1]% of the moves were accepted, reduce delta (only for that parameter)
          delta[foo < thr[1]] <- 0.9*delta[foo < thr[1]]
          # if more than thr[2]% of the moves were accepted, increase delta (only for that parameter)
          delta[foo > thr[2]] <- 1.1*delta[foo > thr[2]]
        }
      }
    }
  }
  return(list(itns=itns, counts=counts, delta=delta))
}
