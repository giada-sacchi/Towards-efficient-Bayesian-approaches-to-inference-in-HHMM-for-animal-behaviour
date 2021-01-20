#### Packages & Scripts ####
if(!require(boot)) install.packages("boot")
library(boot)
if(!require(Rcpp)) install.packages("Rcpp")
library(Rcpp)
if(!require(RcppArmadillo)) install.packages("RcppArmadillo")
library(RcppArmadillo)   
if(!require(invgamma)) install.packages("invgamma")
library(invgamma)
if(!require(numbers)) install.packages("numbers")
library(numbers) # mod
sourceCpp("13253_2017_282_MOESM4_ESM.cpp")


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
            num <- newlikhood + log(dbeta(par.list[[j]][i], alpha[i], beta[i]))
            den <- likhood + log(dbeta(oldpar, alpha[i], beta[i]))
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
            num <- newlikhood + log(dlnorm(par.list[[j]][i], mu[i], sd[i]))
            den <- likhood + log(dlnorm(oldpar, mu[i], sd[i]))
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
            num <- newlikhood + log(dinvgamma(par.list[[j]][i], shape[i], rate[i]))
            den <- likhood + log(dinvgamma(oldpar, shape[i], rate[i]))
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
# Function to compute the sum of logprior distributions of a parameter vector

prior <- function(par.list, mu, sigma, alpha, beta, shape, rate){
  prior <- 0
  for (j in 5:11){ # positions corresponding to target variable parameters
    for (i in 1:3) {
      if (j==11) { 
        prior <- prior +log(dbeta(par.list[[j]][i], alpha[i], beta[i]))
        } else if (is.element(j, c(5,7,9))){
          prior <- prior + log(dlnorm(par.list[[j]][i], mu[i], sigma[i]))
        } else {  
          prior <- prior + log(dinvgamma(par.list[[j]][i], shape[i], rate[i]))
        }
    }
  }
  prior
}
