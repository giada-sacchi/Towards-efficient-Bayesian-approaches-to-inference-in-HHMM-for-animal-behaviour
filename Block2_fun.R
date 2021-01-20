#### Packages & Scripts ####
if(!require(boot)) install.packages("boot")
library(boot)
if(!require(Rcpp)) install.packages("Rcpp")
library(Rcpp)
if(!require(RcppArmadillo)) install.packages("RcppArmadillo")
library(RcppArmadillo)   
if(!require(invgamma)) install.packages("invgamma")
library(invgamma)   # ingamma
if(!require(numbers)) install.packages("numbers")
library(numbers) # mod
sourceCpp("13253_2017_282_MOESM4_ESM.cpp")


##### UPDATEPAR.BLOCK #####
# Function for updating the parameter vector through block updating (Version 2)
updatepar.block2 <- function(par.list, npar, data, n.iter, tuning=F, thr=c(0.2,0.4),
                             likhood,       # current likelihood
                             mu, sigma, alpha, beta, shape, rate, n.tune=n.iter,
                             delta){        # matrix with sd of random walk proposal distribution 
  # Define an array to store sample from posterior distribution for each chain
  itns <- array(0, dim=c(npar, n.iter))
  # counts of the accepted moves for each parameter (rows=parameters, cols=prod states)
  counts <- array(0, dim=c(11,3, n.iter))
  for (t in 1:n.iter){ ## MH iterations
      # Keep a record of the current parameter value being updated 
      oldpar <- par.list[[11]]
      # Propose a new candidate value using a random walk from normal proposal 
      par.list[[11]] <- rnorm(3, par.list[[11]], delta[11,])
      if (all(par.list[[11]]>=0) & all(par.list[[11]]<=1)){ # dw.pi is a probability
        # Compute the log-likelihood (of the move)
        newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
        if(is.nan(newlikhood)){
          A <- 0
        } else {
          # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
          num <- newlikhood + sum(dbeta(par.list[[11]], alpha, beta,log=TRUE))
          den <- likhood + sum(dbeta(oldpar, alpha, beta,log=TRUE))
          A <- min(1, exp(num-den))
        }
      }  else { # if parameter is outside range set the acceptance probability A to 0
        A <- 0 
      }
      if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
        # Accept the move with probability A and store its likelihood
        likhood <- newlikhood
        counts[11,,t] <- 1
      } else {     
        # Reject move and return parameter value to previous value
        par.list[[11]] <- oldpar
      }
    for (j in c(5,7,9)){ # parameter is a mean
        # Keep a record of the current parameter value being updated 
        oldpar <- par.list[[j]]
        # Propose a new candidate value using a random walk from normal proposal 
        par.list[[j]] <- rnorm(3, par.list[[j]], delta[j,])
        if (par.list[[j]][1]>=0 & par.list[[j]][2]>=0 & par.list[[j]][3]>=0){
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
          # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
          num <- newlikhood + sum(dlnorm(par.list[[j]], mu, sd,log=TRUE))
          den <- likhood + sum(dlnorm(oldpar, mu, sd,log=TRUE))
          A <- min(1, exp(num-den))
          } 
        }  else { # Otherwise set the acceptance probability A to 0
          A <- 0 
        }
        if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
          # Accept the move with probability A and store its likelihood
          likhood <- newlikhood
          counts[j,,t] <- 1
        } else {     
          # Reject move and return parameter value to previous value
          par.list[[j]] <- oldpar
        }
    }
    for (j in c(6,8,10)) { # parameter is a sd
        # Keep a record of the current parameter value being updated 
        oldpar <- par.list[[j]]
        # Propose a new candidate value using a random walk from normal proposal 
        par.list[[j]] <- rnorm(3, par.list[[j]], delta[j,])
        if (par.list[[j]][1]>=0 & par.list[[j]][2]>=0 & par.list[[j]][3]>=0){
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
            # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
            num <- newlikhood + sum(dinvgamma(par.list[[j]], shape, rate,log=TRUE))
            den <- likhood + sum(dinvgamma(oldpar, shape, rate,log=TRUE))
            A <- min(1, exp(num-den))
          }
        }  else { # Otherwise set the acceptance probability A to 0
          A <- 0 
        }
        if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
          # Accept the move with probability A and store its likelihood
          likhood <- newlikhood
          counts[j,,t] <- 1
        } else {     
          # Reject move and return parameter value to previous value
          par.list[[j]] <- oldpar
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
  } ## MH iterations
  return(list(itns=itns, counts=counts, delta=delta))
}

